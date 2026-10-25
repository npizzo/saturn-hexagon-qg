const PROJECT_ROOT = dirname(@__DIR__)
environment_digest() = bytes2hex(sha256(read(joinpath(PROJECT_ROOT, "Manifest.toml"))))

function source_digest()
    names = sort(filter(f -> endswith(f, ".jl"), readdir(@__DIR__)))
    return bytes2hex(sha256(join([f*"\n"*read(joinpath(@__DIR__, f), String) for f in names], "\n")))
end

function save_checkpoint(path, sim)
    mkpath(dirname(abspath(path)))
    tmp = path*".partial"
    jldopen(tmp, "w") do f
        f["schema_version"] = 1
        f["source_sha256"] = source_digest()
        f["manifest_sha256"] = environment_digest()
        f["julia_version"] = string(VERSION)
        f["config"] = sim.config
        f["sol"] = Array(sim.prob.sol)
        f["time_s"] = sim.prob.clock.t
        f["step"] = sim.prob.clock.step
        f["convection"] = sim.convection
        f["device"] = sim.device
    end
    mv(tmp, path; force=true)
    return path
end

function compatible_config(c)
    z = deepcopy(c)
    delete!(z["numerics"], "duration_s")
    return z
end

function restore_checkpoint!(sim, path)
    d = JLD2.load(path)
    d["schema_version"] == 1 || error("Unsupported checkpoint schema")
    d["source_sha256"] == source_digest() || error("Checkpoint source differs; restart requires the original source")
    d["manifest_sha256"] == environment_digest() || error("Checkpoint dependency Manifest differs; restore the original environment")
    d["julia_version"] == string(VERSION) || error("Restart requires identical Julia version for RNG/serialization reproducibility")
    compatible_config(d["config"]) == compatible_config(sim.config) || error("Restart config differs (only duration_s may change)")
    size(d["sol"]) == size(sim.prob.sol) || error("Checkpoint grid mismatch")
    d["time_s"] < sim.config["numerics"]["duration_s"] || error("Restart duration must extend beyond checkpoint time")
    copyto!(sim.prob.sol, d["sol"])
    sim.prob.clock.t = d["time_s"]; sim.prob.clock.step = d["step"]
    old = d["convection"]
    sim.convection.next_event = old.next_event
    sim.convection.events = old.events
    sim.convection.rng = old.rng
    sim.convection.cfg = old.cfg
    QG.updatevars!(sim.prob)
    return sim
end

"Physical regime indicators inside the untapered science disk."
function physical_diagnostics(sim, f=fields(sim))
    c = sim.config
    ld = c["physics"]["deformation_radius_m"]
    kap2 = isinf(ld) ? 0.0 : inv(ld^2)
    f0 = 2c["planet"]["rotation_rate_s"]
    gamma = f0/c["planet"]["radius_m"]^2
    r2 = [xx^2+yy^2 for xx in sim.x, yy in sim.y]
    disk = r2 .<= c["boundary"]["inner_radius_m"]^2
    # s is Helmholtz PV, not relative vorticity. The constant f0 omitted
    # from prognostic PV still belongs in this local Coriolis denominator.
    zeta = f.s .+ kap2 .* f.psi
    local_f = f0 .- gamma/2 .* r2
    maxro = maximum(abs.(zeta[disk] ./ local_f[disk]))
    if isinf(ld)
        return (; max_abs_hprime_over_H=NaN, min_h_over_H=NaN, max_abs_zeta_over_f=maxro)
    end
    chi = kap2/f0 .* (f.psi[disk] .- sim.psid[disk])
    return (; max_abs_hprime_over_H=maximum(abs, chi), min_h_over_H=1+minimum(chi), max_abs_zeta_over_f=maxro)
end

function diagnostics(sim)
    f = fields(sim); c = sim.config; o = c["output"]
    pol = polar_diagnostics(f.u, f.v, f.totalpv, f.psi, f.s, sim.x, sim.y;
        jet_radius_m=c["initial"]["jet_radius_m"], jet_halfwidth_m=o["jet_search_halfwidth_m"],
        m_max=o["m_max"], ntheta=o["ntheta"], nr=o["nr"])
    inv = invariants(sim)
    physical = physical_diagnostics(sim, f)
    header = ["time_s", "step", "energy_m2_s2", "total_pv_enstrophy_s2", "mean_s_s", "mean_jet_radius_m", "mean_peak_speed_m_s", "jet_coverage"]
    vals = Any[sim.prob.clock.t, sim.prob.clock.step, inv.energy, inv.total_enstrophy, inv.mean_s, pol.meanjet, pol.mean_peak_speed, pol.coverage]
    append!(header, string.(collect(keys(physical))))
    append!(vals, collect(values(physical)))
    for m in 1:o["m_max"]
        append!(header, ["m$(m)_amplitude_m", "m$(m)_real_m", "m$(m)_imag_m"])
        z = pol.mode_coefficients[m]
        append!(vals, [pol.mode_amplitudes[m], real(z), imag(z)])
    end
    return (; header, vals, polar=pol)
end

function append_csv(path, header, values)
    fresh = !isfile(path)
    open(path, "a") do io
        fresh && println(io, join(header, ','))
        println(io, join(values, ','))
    end
end

function write_profile(path, pol)
    open(path, "w") do io
        println(io, "radius_m,u_theta_m_s")
        for i in eachindex(pol.radii)
            println(io, join((pol.radii[i], pol.mean_utheta[i]), ','))
        end
    end
end

function write_radial_diagnostics(path, pol)
    open(path, "w") do io
        println(io, "radius_m,u_theta_m_s,u_radial_m_s,total_pv_s,eddy_ur_utheta_m2_s2")
        for i in eachindex(pol.radii)
            println(io, join((pol.radii[i], pol.mean_utheta[i], pol.mean_ur[i], pol.mean_pv[i], pol.eddy_uv[i]), ','))
        end
    end
end

function write_snapshot(path, sim)
    f = fields(sim)
    jldopen(path, "w") do out
        out["time_s"] = sim.prob.clock.t
        out["step"] = sim.prob.clock.step
        for k in keys(f)
            out[string(k)] = getfield(f, k)
        end
    end
end

"Run into a NEW directory; never overwrite an earlier experiment."
function run_experiment(config_path; output, device="cpu", resume=nothing)
    c = load_config(config_path)
    ispath(output) && error("Output already exists: $output; choose a new directory (also for resume)")
    sim = build_simulation(c; device)
    resume !== nothing && restore_checkpoint!(sim, resume)
    mkpath(output)
    cp(config_path, joinpath(output, "config.toml"))
    for name in ("Project.toml", "Manifest.toml")
        src = joinpath(PROJECT_ROOT, name)
        isfile(src) && cp(src, joinpath(output, name))
    end
    cp(@__DIR__, joinpath(output, "src"))
    metadata = Dict("schema_version"=>1, "created_utc"=>string(now(UTC)), "julia_version"=>string(VERSION),
        "source_sha256"=>source_digest(), "manifest_sha256"=>environment_digest(), "device"=>device, "threads"=>Threads.nthreads(),
        "axes"=>"arrays indexed [x,y]; x=r*cos(theta), y=r*sin(theta); u_theta=-u*sin(theta)+v*cos(theta)",
        "units"=>"SI; planar distance in m; time s; velocity m/s; psi m2/s; s and totalpv 1/s",
        "pv_definition"=>"s=(Laplacian-LD^-2)*psi; totalpv=s+B; B=LD^-2*psid+continued gamma-plane PV",
        "deep_profile"=>get(c["physics"], "deep_profile", "annular"),
        "physical_diagnostics"=>"r<=boundary.inner_radius_m; hprime/H=(psi-psid)/(2*rotation_rate_s*LD^2), undefined at infinite LD; h/H=1+hprime/H; zeta=s+psi/LD^2; f(r)=2*rotation_rate_s-(2*rotation_rate_s/radius_m^2)*r^2/2",
        "resume_from"=>resume === nothing ? "" : abspath(resume), "hyperviscosity"=>sim.nu)
    open(joinpath(output, "metadata.toml"), "w") do io; TOML.print(io, metadata); end
    jldsave(joinpath(output, "domain.jld2"); x=sim.x, y=sim.y, background_pv=sim.background, psi_deep=sim.psid, reference_s=sim.reference_s)
    o = c["output"]; t0 = sim.prob.clock.t
    nextdiag = t0; nextsnap = t0; nextcheck = t0+o["checkpoint_interval_s"]
    lastdiag = -Inf
    function record!()
        d = diagnostics(sim)
        append_csv(joinpath(output, "diagnostics.csv"), d.header, d.vals)
        write_profile(joinpath(output, "profile_latest.csv"), d.polar)
        write_radial_diagnostics(joinpath(output, "radial_diagnostics_latest.csv"), d.polar)
        return nothing
    end
    record!(); lastdiag = t0; nextdiag += o["diagnostic_interval_s"]
    write_snapshot(joinpath(output, @sprintf("snapshot_%09d.jld2", sim.prob.clock.step)), sim)
    nextsnap += o["snapshot_interval_s"]
    try
        while sim.prob.clock.t < c["numerics"]["duration_s"] - 8eps(c["numerics"]["duration_s"])
            event = advance!(sim)
            if event !== nothing
                append_csv(joinpath(output, "events.csv"), string.(collect(keys(event))), collect(values(event)))
            end
            t = sim.prob.clock.t
            if t >= nextdiag
                record!(); lastdiag = t
                nextdiag = (floor((t-t0)/o["diagnostic_interval_s"])+1)*o["diagnostic_interval_s"]+t0
            end
            if t >= nextsnap
                write_snapshot(joinpath(output, @sprintf("snapshot_%09d.jld2", sim.prob.clock.step)), sim)
                nextsnap = (floor((t-t0)/o["snapshot_interval_s"])+1)*o["snapshot_interval_s"]+t0
            end
            if t >= nextcheck
                save_checkpoint(joinpath(output, "checkpoint.jld2"), sim)
                nextcheck = (floor((t-t0)/o["checkpoint_interval_s"])+1)*o["checkpoint_interval_s"]+t0
            end
        end
        sim.prob.clock.t > lastdiag && record!()
        write_snapshot(joinpath(output, @sprintf("snapshot_%09d.jld2", sim.prob.clock.step)), sim)
        save_checkpoint(joinpath(output, "checkpoint.jld2"), sim)
        open(joinpath(output, "COMPLETED"), "w") do io; println(io, sim.prob.clock.t); end
    catch err
        open(joinpath(output, "FAILED.txt"), "w") do io; showerror(io, err); end
        rethrow()
    end
    return sim
end
