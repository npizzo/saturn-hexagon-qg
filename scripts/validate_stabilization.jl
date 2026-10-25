module StabilizationValidation

using SaturnPolarQG, CUDA, LinearAlgebra, Random, SHA, TOML, Dates, Printf, JLD2

const FIELD_RTOL = 1e-8
const BASE_THICKNESS_ATOL = 1e-11
const FIELD_KEYS = (:s, :psi, :u, :v, :totalpv)

function options(args)
    length(args) == 3 && args[2] == "--output" ||
        error("Usage: scripts/validate_stabilization.jl CONFIG --output NEW_DIRECTORY")
    return abspath(args[1]), abspath(args[3])
end

function field_errors(a, b; tolerance=FIELD_RTOL)
    errors = Dict{String,Any}()
    for key in FIELD_KEYS
        x, y = getproperty(a,key), getproperty(b,key)
        size(x) == size(y) || error("CPU/GPU $key dimensions differ")
        all(isfinite,x) && all(isfinite,y) || error("Nonfinite CPU/GPU $key")
        relative = norm(x-y)/max(norm(x),eps(Float64))
        isfinite(relative) && relative <= tolerance || error("CPU/GPU $key relative L2 difference $relative exceeds $tolerance")
        errors[string(key)] = Dict("relative_l2"=>relative,"maximum_absolute"=>maximum(abs,x-y))
    end
    return errors
end

forcing_state(sim) = (; time_s=sim.prob.clock.t, step=sim.prob.clock.step,
    events=sim.convection.events, next_event_s=sim.convection.next_event,
    rng_probe=rand(deepcopy(sim.convection.rng),UInt64,16))

function check_forcing_states(a, b, enabled)
    a == b || error("CPU/GPU clock, step, event schedule, or random stream differs")
    if enabled
        a.events > 0 || error("Enabled-convection validation did not exercise an event")
    else
        a.events == 0 && a.next_event_s == Inf || error("Disabled convection must have zero events and an infinite next-event time")
    end
    return nothing
end

function snapshot_fields(path)
    d = JLD2.load(path)
    return NamedTuple{FIELD_KEYS}(Tuple(d[string(k)] for k in FIELD_KEYS))
end

function full_match_applicable(c)
    d = c["physics"]
    profile = get(d,"deep_profile","annular")
    matched = profile == "matched_upper" ? get(d,"deep_match_fraction",1.0) == 1 :
        profile == "matched_components" ? get(d,"deep_jet_fraction",NaN) == 1 && get(d,"deep_cyclone_fraction",NaN) == 1 : false
    return matched && isfinite(d["deformation_radius_m"])
end

function check_full_match_base(c)
    if !full_match_applicable(c)
        return Dict{String,Any}("applicable"=>false,
            "reason"=>"Zero base thickness is asserted only for finite LD with every upper-flow component fully matched.")
    end
    baseconfig = deepcopy(c)
    baseconfig["initial"]["perturbation"] = "none"
    baseconfig["convection"]["enabled"] = false
    checks = Dict{String,Any}("applicable"=>true,"absolute_thickness_tolerance"=>BASE_THICKNESS_ATOL)
    for device in ("cpu","gpu")
        sim = build_simulation(baseconfig; device)
        p = SaturnPolarQG.physical_diagnostics(sim)
        isfinite(p.max_abs_hprime_over_H) && p.max_abs_hprime_over_H < BASE_THICKNESS_ATOL ||
            error("$device unperturbed full-match thickness cancellation failed")
        abs(p.min_h_over_H-1) < BASE_THICKNESS_ATOL || error("$device full-match base thickness differs from H")
        checks[device] = Dict(string(k)=>v for (k,v) in pairs(p))
    end
    return checks
end

function mark_failed(output, err)
    if isdir(output)
        # An integration finishing is not sufficient for a validation pass.
        rm(joinpath(output,"COMPLETED"); force=true)
        rm(joinpath(output,"CPU_GPU_VALIDATION_PASSED"); force=true)
        open(joinpath(output,"FAILED.txt"),"w") do io
            showerror(io,err)
        end
    end
    return nothing
end

function validate(config_path, output)
    c = load_config(config_path)
    c["experiment"]["purpose"] == "smoke" || error("Use a purpose=smoke validation configuration")
    c["numerics"]["duration_s"] <= 20_000 || error("Validation duration_s must be <=20000")
    ispath(output) && error("Output already exists: $output; choose a new directory")
    CUDA.functional(true) || error("Stabilization CPU/GPU validation requires working CUDA; no CPU fallback")
    CUDA.allowscalar(false)
    CUDA.versioninfo()
    report = Dict{String,Any}(
        "schema_version"=>1,"created_utc"=>string(now(UTC)),"julia_version"=>string(VERSION),
        "source_sha256"=>SaturnPolarQG.source_digest(),"manifest_sha256"=>SaturnPolarQG.environment_digest(),
        "config_sha256"=>bytes2hex(sha256(read(config_path))),
        "validator_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "relative_l2_tolerance"=>FIELD_RTOL,"convection_enabled"=>c["convection"]["enabled"],
        "deep_profile"=>get(c["physics"],"deep_profile","annular"),
        "initial_profile"=>get(c["initial"],"profile","annular"),
        "scope"=>"Short-time CPU/GPU agreement and applicable full-match base cancellation; not physical stability, spatial convergence, or long-time chaotic trajectory equivalence.")
    try
        cpu = build_simulation(c; device="cpu")
        initial_cpu = fields(cpu)
        # Retain the ordinary run layout for existing collectors and galleries.
        gpu = run_experiment(config_path; output,device="gpu")
        completion_text = read(joinpath(output,"COMPLETED"),String)
        rm(joinpath(output,"COMPLETED")) # publish completion only after validation
        cp(@__FILE__,joinpath(output,"validate_stabilization.jl"))
        initial_gpu = snapshot_fields(joinpath(output,"snapshot_000000000.jld2"))
        report["initial_errors"] = field_errors(initial_cpu,initial_gpu)
        while cpu.prob.clock.t < c["numerics"]["duration_s"]
            advance!(cpu)
        end
        report["final_errors"] = field_errors(fields(cpu),fields(gpu))
        a, b = forcing_state(cpu), forcing_state(gpu)
        check_forcing_states(a,b,c["convection"]["enabled"])
        abs(invariants(gpu).mean_s) < 1e-15 || error("GPU mean-PV compensation check failed")
        report["final_clock"] = Dict("time_s"=>b.time_s,"step"=>b.step,"events"=>b.events,"next_event_s"=>b.next_event_s)
        report["exact_checks"] = ["final time","step count","event count","next-event time","16-draw UInt64 RNG probe"]
        report["full_match_base"] = check_full_match_base(c)
        report["final_gpu_physical_indicators"] = Dict(string(k)=>v for (k,v) in pairs(SaturnPolarQG.physical_diagnostics(gpu)))
        report["passed"] = true
        open(joinpath(output,"cpu_gpu_validation.toml.partial"),"w") do io
            TOML.print(io,report; sorted=true)
        end
        mv(joinpath(output,"cpu_gpu_validation.toml.partial"),joinpath(output,"cpu_gpu_validation.toml"))
        write(joinpath(output,"CPU_GPU_VALIDATION_PASSED"),"Short-time CPU/GPU validation passed.\n")
        write(joinpath(output,"COMPLETED"),completion_text)
        for key in sort(collect(keys(report["final_errors"])))
            @printf("%s final CPU/GPU relative L2 difference: %.6e\n",key,report["final_errors"][key]["relative_l2"])
        end
        println("STABILIZATION CPU/GPU VALIDATION PASSED; results: ",output)
        return report
    catch err
        mark_failed(output,err)
        rethrow()
    end
end

function main(args)
    if args == ["--help"] || args == ["-h"]
        println("Usage: scripts/validate_stabilization.jl CONFIG --output NEW_DIRECTORY\nRequires CUDA; compares short CPU/GPU trajectories and preserves ordinary GPU run artifacts.")
        return
    end
    config,output = options(args)
    validate(config,output)
end

end # module

abspath(PROGRAM_FILE) == (@__FILE__) && StabilizationValidation.main(ARGS)
