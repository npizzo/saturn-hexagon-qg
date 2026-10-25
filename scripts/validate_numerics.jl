using SaturnPolarQG
using SaturnPolarQG.PolarDiagnostics: radial_anomaly, local_pv_anomaly
using Statistics
using Printf

const VALIDATION_PROJECT = dirname(@__DIR__)

function validation_config(nx, dt)
    c = load_config(joinpath(VALIDATION_PROJECT, "configs", "smoke.toml"))
    c["grid"]["nx"] = nx
    c["convection"]["enabled"] = false
    c["boundary"]["sponge_rate_s"] = 0.0
    c["numerics"]["cutoff_damping_time_s"] = Inf
    c["numerics"]["dt_s"] = dt
    c["numerics"]["duration_s"] = 4000.0
    return c
end

function axisymmetric_limiter_residual(c)
    axisymmetric = deepcopy(c)
    axisymmetric["initial"]["perturbation"] = "none"
    sim = build_simulation(axisymmetric; device="cpu")
    f = fields(sim)
    anomaly = radial_anomaly(f.totalpv, sim.x, sim.y)
    radius = [hypot(x,y) for x in sim.x, y in sim.y]
    a, v = c["initial"], c["convection"]
    jet = (radius .>= a["jet_radius_m"]-a["jet_width_m"]) .&
          (radius .<= a["jet_radius_m"]+a["jet_width_m"])
    convection = (radius .>= v["r_min_m"]) .& (radius .<= v["r_max_m"])
    region_stats(mask) = (sqrt(mean(abs2, anomaly[mask])), maximum(abs, anomaly[mask]))
    jet_rms, jet_max = region_stats(jet)
    convection_rms, convection_max = region_stats(convection)
    qmax = v["q_max_s"]
    # Fixed, equally spaced r² and theta samples approximate an area-weighted
    # annular statistic. Compare both limiters at exactly these same points.
    sample_nr, sample_ntheta = 16, 48
    dx, dy = sim.x[2]-sim.x[1], sim.y[2]-sim.y[1]
    function sampled_stats(lo,hi)
        old, current = Float64[], Float64[]
        for r2 in range(lo^2,hi^2; length=sample_nr), k in 0:sample_ntheta-1
            theta = 2pi*k/sample_ntheta
            xp, yp = sqrt(r2)*cos(theta), sqrt(r2)*sin(theta)
            push!(old,SaturnPolarQG.Convection._sample(anomaly,xp,yp,sim.x,sim.y,dx,dy))
            push!(current,local_pv_anomaly(f.totalpv,sim.x,sim.y,xp,yp))
        end
        return (; legacy_rms=sqrt(mean(abs2,old)), legacy_max=maximum(abs,old),
            local_rms=sqrt(mean(abs2,current)), local_max=maximum(abs,current))
    end
    sampled_jet = sampled_stats(a["jet_radius_m"]-a["jet_width_m"],a["jet_radius_m"]+a["jet_width_m"])
    sampled_convection = sampled_stats(v["r_min_m"],v["r_max_m"])
    return (; radial_bins=max(16, length(sim.x)÷2), qmax_s=qmax,
        legacy_grid_jet_rms_s=jet_rms, legacy_grid_jet_max_s=jet_max,
        legacy_grid_jet_rms_over_qmax=jet_rms/qmax, legacy_grid_jet_max_over_qmax=jet_max/qmax,
        legacy_grid_convection_rms_s=convection_rms, legacy_grid_convection_max_s=convection_max,
        legacy_grid_convection_rms_over_qmax=convection_rms/qmax,
        legacy_grid_convection_max_over_qmax=convection_max/qmax,
        sample_nr,sample_ntheta,
        legacy_sampled_jet_rms_over_qmax=sampled_jet.legacy_rms/qmax,
        legacy_sampled_jet_max_over_qmax=sampled_jet.legacy_max/qmax,
        local_sampled_jet_rms_over_qmax=sampled_jet.local_rms/qmax,
        local_sampled_jet_max_over_qmax=sampled_jet.local_max/qmax,
        legacy_sampled_convection_rms_over_qmax=sampled_convection.legacy_rms/qmax,
        legacy_sampled_convection_max_over_qmax=sampled_convection.legacy_max/qmax,
        local_sampled_convection_rms_over_qmax=sampled_convection.local_rms/qmax,
        local_sampled_convection_max_over_qmax=sampled_convection.local_max/qmax)
end

function conservation_check(c)
    sim = build_simulation(c; device="cpu")
    initial = invariants(sim)
    initial_polar = SaturnPolarQG.diagnostics(sim).polar
    max_energy_drift, max_enstrophy_drift, max_mean_drift = 0.0, 0.0, 0.0
    final = initial
    while sim.prob.clock.t < c["numerics"]["duration_s"]
        advance!(sim)
        final = invariants(sim)
        all(isfinite, (final.energy, final.total_enstrophy, final.mean_s)) || error("Nonfinite invariant")
        max_energy_drift = max(max_energy_drift, abs(final.energy/initial.energy-1))
        max_enstrophy_drift = max(max_enstrophy_drift, abs(final.total_enstrophy/initial.total_enstrophy-1))
        max_mean_drift = max(max_mean_drift, abs(final.mean_s-initial.mean_s))
    end
    return (; nx=c["grid"]["nx"], dx_m=c["grid"]["length_m"]/c["grid"]["nx"],
        dt_s=c["numerics"]["dt_s"], duration_s=sim.prob.clock.t, steps=sim.prob.clock.step,
        energy_initial_m2_s2=initial.energy, energy_final_m2_s2=final.energy,
        energy_relative_final=final.energy/initial.energy-1,
        energy_max_abs_relative=max_energy_drift,
        total_enstrophy_initial_s2=initial.total_enstrophy,
        total_enstrophy_final_s2=final.total_enstrophy,
        total_enstrophy_relative_final=final.total_enstrophy/initial.total_enstrophy-1,
        total_enstrophy_max_abs_relative=max_enstrophy_drift,
        mean_s_max_abs_change_s=max_mean_drift,
        initial_jet_coverage=initial_polar.coverage,
        initial_m2_amplitude_m=initial_polar.mode_amplitudes[2],
        initial_m4_amplitude_m=initial_polar.mode_amplitudes[4],
        initial_m6_amplitude_m=initial_polar.mode_amplitudes[6],
        initial_m10_amplitude_m=initial_polar.mode_amplitudes[10],
        initial_m12_amplitude_m=initial_polar.mode_amplitudes[12])
end

function validation_main(args)
    length(args) in (0,2) || error("Usage: julia --project=. scripts/validate_numerics.jl [--output NEW_CSV]")
    !isempty(args) && args[1] != "--output" && error("Only --output is supported")
    path = isempty(args) ? joinpath(VALIDATION_PROJECT,"validation","numerical_checks.csv") : abspath(args[2])
    ispath(path) && error("Refusing to overwrite $path; supply a new --output path")
    mkpath(dirname(path))
    source_hash = SaturnPolarQG.source_digest()
    manifest_hash = SaturnPolarQG.environment_digest()
    println("CPU numerical validation: unforced, no sponge or viscosity, fixed physical profiles.")
    println("Each resolution advances the m=6 perturbed state for 4000 s; limiter residuals use a separate axisymmetric state at t=0.")
    println("The 64² case is also repeated with dt=10 s to check timestep sensitivity.")
    # Stream rows as cases finish; .partial is never represented as a complete result.
    partial = path*".partial"
    residuals = Dict{Int,NamedTuple}()
    open(partial,"w") do io
        for (index,(nx,dt)) in enumerate(((64,20.0),(128,20.0),(256,20.0),(64,10.0)))
            c = validation_config(nx,dt)
            residual = get!(residuals,nx) do
                axisymmetric_limiter_residual(c)
            end
            conservation = conservation_check(c)
            row = merge(conservation, residual,
                (; julia_version=string(VERSION), julia_threads=Threads.nthreads(),
                   source_sha256=source_hash, manifest_sha256=manifest_hash))
            index == 1 && println(io,join(string.(keys(row)),','))
            println(io,join(values(row),',')); flush(io)
            @printf("nx=%d, dt=%.0f s: max |ΔE/E|=%.3e, max |ΔZ/Z|=%.3e; initial axisymmetric LOCAL q′ RMS/qmax: jet %.3e, convection %.3e\n",
                nx,dt,row.energy_max_abs_relative,row.total_enstrophy_max_abs_relative,
                row.local_sampled_jet_rms_over_qmax,row.local_sampled_convection_rms_over_qmax)
        end
    end
    SaturnPolarQG.source_digest() == source_hash || error("Source changed during validation; results left as .partial")
    SaturnPolarQG.environment_digest() == manifest_hash || error("Manifest changed during validation; results left as .partial")
    mv(partial,path)
    println("Saved ",path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    validation_main(ARGS)
end
