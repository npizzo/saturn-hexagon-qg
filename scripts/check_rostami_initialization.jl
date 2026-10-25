#!/usr/bin/env julia
# Audit the actual projected circular state; does not test nonlinear stability.
using SaturnPolarQG, FFTW, LinearAlgebra, Statistics, TOML, SHA

function audit_initialization(config_path; steps=3)
    c=load_config(config_path)
    get(c["initial"],"profile","annular")=="rostami2017" || error("A Rostami profile is required")
    original_perturbation=c["initial"]["perturbation"]
    c["initial"]["perturbation"]="none"
    c["convection"]["enabled"]=false
    sim=build_simulation(c;device="cpu")
    f=fields(sim); g=sim.prob.grid; a=c["initial"]
    radius=[hypot(x,y) for x in sim.x,y in sim.y]
    disk=radius .<= c["boundary"]["inner_radius_m"]
    jet=abs.(radius.-a["jet_radius_m"]) .<= 2a["jet_width_m"]
    core=radius .<= 3a["cyclone_radius_m"]
    qh=rfft(f.totalpv)
    qx=irfft(im.*g.kr.*qh,length(sim.x)); qy=irfft(im.*g.l.*qh,length(sim.x))
    residual=f.u.*qx.+f.v.*qy
    prof=SaturnPolarQG.rostami_profile(a)
    reference_u=zeros(size(radius));reference_v=similar(reference_u)
    for j in eachindex(sim.y),i in eachindex(sim.x)
        r=radius[i,j]
        speed=SaturnPolarQG.rostami_velocity(r,prof)
        reference_u[i,j]=r==0 ? 0.0 : -speed*sim.y[j]/r
        reference_v[i,j]=r==0 ? 0.0 : speed*sim.x[i]/r
    end
    tref=2pi*a["jet_radius_m"]/abs(a["jet_speed_m_s"])
    qref=abs(a["jet_speed_m_s"])/a["jet_width_m"]
    zones=Dict{String,Any}()
    for (name,mask) in (("science_disk",disk),("jet",jet),("core",core))
        rms=sqrt(mean(abs2,residual[mask]))
        velocity_error=sqrt(mean((f.u[mask].-reference_u[mask]).^2 .+(f.v[mask].-reference_v[mask]).^2))
        zones[name]=Dict("rms_inviscid_residual_s2"=>rms,
            "maximum_inviscid_residual_s2"=>maximum(abs,residual[mask]),
            "dimensionless_residual"=>tref*rms/qref,
            "rms_velocity_projection_error_m_s"=>velocity_error)
    end
    # Explicit finite-step checks exercise the production CFL/background/damping guards.
    # Hyperviscosity remains enabled: this trajectory is not an inviscid equilibrium test.
    for _ in 1:steps
        advance!(sim)
    end
    after=fields(sim)
    all(isfinite,after.s) || error("Nonfinite evolved PV")
    change=norm(after.psi-f.psi)/max(norm(f.psi),eps())
    physical=SaturnPolarQG.physical_diagnostics(sim,f)
    return Dict("config"=>basename(config_path),"config_sha256"=>bytes2hex(sha256(read(config_path))),
        "source_sha256"=>SaturnPolarQG.source_digest(),"manifest_sha256"=>SaturnPolarQG.environment_digest(),
        "scope"=>"Circularized actual-grid CPU initialization audit. Not a GPU comparison, a long integration, or a stability certificate.",
        "original_perturbation"=>original_perturbation,"audited_perturbation"=>"none",
        "nx"=>c["grid"]["nx"],"core_regularization_m"=>a["core_regularization_m"],
        "jet_circulation_time_s"=>tref,"pv_scale_s"=>qref,"zones"=>zones,
        "initial_physical_indicators"=>Dict(string(k)=>v for (k,v) in pairs(physical)),
        "advanced_steps"=>steps,"advanced_time_s"=>sim.prob.clock.t,"short_psi_relative_change"=>change,
        "interpretation"=>"Compare residuals across grids at fixed core scale. Core smoothing changes the physical profile. QG validity requires small Rossby number and interface displacement, not just small residuals.")
end

function main(args)
    length(args)>=3 && args[1]=="--output" || error("Usage: check_rostami_initialization.jl --output NEW_DIRECTORY CONFIG...")
    out=abspath(args[2]); ispath(out) && error("Choose a new output directory")
    configs=abspath.(args[3:end]); all(isfile,configs)||error("Missing config")
    allunique(basename.(configs))||error("Config basenames must be unique")
    mkpath(out)
    for path in configs
        report=audit_initialization(path)
        open(joinpath(out,splitext(basename(path))[1]*".toml"),"w") do io
            TOML.print(io,report;sorted=true)
        end
        println(basename(path),": N=",report["nx"],"; disk residual=",report["zones"]["science_disk"]["dimensionless_residual"],
            "; max |zeta/f|=",report["initial_physical_indicators"]["max_abs_zeta_over_f"])
        GC.gc()
    end
end

abspath(PROGRAM_FILE)==(@__FILE__) && main(ARGS)
