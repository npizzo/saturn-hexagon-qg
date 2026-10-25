using SaturnPolarQG, CUDA, LinearAlgebra, Random, SHA, TOML, Dates, Printf, JLD2

const RESTART_RTOL = 1e-12

function restart_options(args)
    isempty(args) && error("Usage: scripts/validate_gpu_restart.jl CONFIG --output NEW_DIRECTORY [--device gpu|cpu]")
    iseven(length(args)-1) || error("Each option requires a value")
    opts = Dict{String,String}()
    for i in 2:2:length(args)
        args[i] in ("--output", "--device") || error("Unknown option $(args[i])")
        haskey(opts, args[i]) && error("Duplicate option $(args[i])")
        opts[args[i]] = args[i+1]
    end
    haskey(opts, "--output") || error("--output is required")
    device = get(opts, "--device", "gpu")
    device in ("cpu", "gpu") || error("device must be cpu or gpu")
    return abspath(args[1]), abspath(opts["--output"]), device
end

# Probe a copy: validation must never consume the simulation's random stream.
rng_probe(sim) = rand(deepcopy(sim.convection.rng), UInt64, 32)
clock_and_forcing(sim) = (; time_s=sim.prob.clock.t, step=sim.prob.clock.step,
    events=sim.convection.events, next_event_s=sim.convection.next_event)

function array_error(a, b)
    all(isfinite, a) && all(isfinite, b) || error("Nonfinite state in restart comparison")
    return Dict("relative_l2"=>norm(a-b)/max(norm(a), eps(Float64)),
                "maximum_absolute"=>maximum(abs, a-b))
end

function advance_record!(sim, clocks, events)
    event = advance!(sim)
    push!(clocks, clock_and_forcing(sim))
    if event !== nothing
        # Positions and event times depend only on the saved random stream.
        push!(events, (; event.event, event.time_s, event.next_event_s,
            event.x_m, event.y_m, event.cyclone_x_m, event.cyclone_y_m,
            event.anticyclone_x_m, event.anticyclone_y_m, event.injected))
    end
    return event
end

function validate_restart(config_path, output, device)
    c = load_config(config_path)
    c["convection"]["enabled"] || error("Restart validation requires enabled convection")
    duration = c["numerics"]["duration_s"]
    duration <= 20_000 || error("Use a short configuration (duration_s <= 20000), such as configs/unity_benchmark.toml")
    ispath(output) && error("Output already exists: $output; choose a new directory")
    if device == "gpu"
        CUDA.functional(true) || error("GPU restart validation requires working CUDA; no CPU fallback")
        CUDA.allowscalar(false)
        CUDA.versioninfo()
    end
    mkpath(dirname(output))
    mkdir(output)
    cp(config_path, joinpath(output, "config.toml"))
    cp(@__FILE__, joinpath(output, "validate_gpu_restart.jl"))
    for name in ("Project.toml", "Manifest.toml", "src")
        cp(joinpath(SaturnPolarQG.PROJECT_ROOT, name), joinpath(output, name))
    end
    report = Dict{String,Any}(
        "schema_version"=>1, "created_utc"=>string(now(UTC)), "device"=>device,
        "julia_version"=>string(VERSION), "threads"=>Threads.nthreads(),
        "source_sha256"=>SaturnPolarQG.source_digest(),
        "manifest_sha256"=>SaturnPolarQG.environment_digest(),
        "config_sha256"=>bytes2hex(sha256(read(config_path))),
        "validator_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "relative_l2_tolerance"=>RESTART_RTOL,
        "scope"=>"Fresh model and serialized checkpoint, same device and allocation; no cross-node or long-time chaotic equivalence claim.")
    checkpoint = joinpath(output, "checkpoint.jld2")
    try
        reference = build_simulation(c; device)
        before_events = NamedTuple[]
        before_clocks = NamedTuple[]
        # Do not set `until` to the split time: that would change RK steps.
        # Save at the first ordinary step/event boundary at or after halfway.
        while reference.prob.clock.t < duration/2
            advance_record!(reference, before_clocks, before_events)
        end
        split = clock_and_forcing(reference)
        split.time_s < duration || error("No interior checkpoint boundary; use a longer short configuration")
        count(e->e.injected, before_events) > 0 || error("No nonzero convection injection before checkpoint")
        split_rng = rng_probe(reference)
        split_cfg = deepcopy(reference.convection.cfg)
        split_sol = Array(reference.prob.sol)
        save_checkpoint(checkpoint, reference)
        JLD2.load(checkpoint, "sol") == split_sol || error("Checkpoint spectral serialization is not exact")
        # updatevars! projects the loaded state. The uninterrupted branch's
        # next check_step does the identical projection. Compare retained modes
        # here without touching the live reference or accepting serialization loss.
        SaturnPolarQG.QG.FourierFlows.dealias!(split_sol, reference.prob.grid)
        @printf("Checkpoint saved at ordinary boundary t=%.12g s, step=%d, events=%d\n", split.time_s, split.step, split.events)
        flush(stdout)

        reference_clocks = NamedTuple[]
        reference_events = NamedTuple[]
        while reference.prob.clock.t < duration
            advance_record!(reference, reference_clocks, reference_events)
        end
        count(e->e.injected, reference_events) > 0 || error("No nonzero convection injection after checkpoint")
        final_state = clock_and_forcing(reference)
        final_rng = rng_probe(reference)
        final_sol = Array(reference.prob.sol)
        final_fields = fields(reference)
        # Release the reference before rebuilding; no integrator work arrays
        # or live RNG object are reused by the checkpoint-loaded trajectory.
        reference = nothing
        GC.gc()

        resumed = restore_checkpoint!(build_simulation(c; device), checkpoint)
        clock_and_forcing(resumed) == split || error("Restored clock/step/event schedule differs")
        rng_probe(resumed) == split_rng || error("Restored random stream differs")
        resumed.convection.cfg == split_cfg || error("Restored convection configuration differs")
        restored_error = array_error(split_sol, Array(resumed.prob.sol))
        restored_error["relative_l2"] <= RESTART_RTOL || error("Restored spectral state differs")
        resumed_clocks = NamedTuple[]
        resumed_events = NamedTuple[]
        while resumed.prob.clock.t < duration
            advance_record!(resumed, resumed_clocks, resumed_events)
        end
        resumed_clocks == reference_clocks || error("Post-restart timestep/event schedule differs")
        resumed_events == reference_events || error("Post-restart random event locations or injection flags differ")
        clock_and_forcing(resumed) == final_state || error("Final clock/step/event schedule differs")
        rng_probe(resumed) == final_rng || error("Final random stream differs")
        errors = Dict("spectral_s"=>array_error(final_sol, Array(resumed.prob.sol)))
        resumed_fields = fields(resumed)
        for key in keys(final_fields)
            errors[string(key)] = array_error(getfield(final_fields, key), getfield(resumed_fields, key))
        end
        for key in sort(collect(keys(errors)))
            @printf("%s restart relative L2 difference: %.6e; maximum absolute: %.6e\n", key,
                errors[key]["relative_l2"], errors[key]["maximum_absolute"])
        end
        maximum(e["relative_l2"] for e in values(errors)) <= RESTART_RTOL || error("Restart field tolerance exceeded")
        report["passed"] = true
        report["checkpoint"] = Dict(string(k)=>v for (k,v) in pairs(split))
        report["final"] = Dict(string(k)=>v for (k,v) in pairs(final_state))
        report["nonzero_injections_before_checkpoint"] = count(e->e.injected, before_events)
        report["nonzero_injections_after_checkpoint"] = count(e->e.injected, reference_events)
        report["restored_spectral_error"] = restored_error
        report["final_errors"] = errors
        report["checkpoint_spectral_comparison"] = "Exact disk serialization, then the library's usual dealias projection before comparing restored retained modes."
        report["exact_checks"] = ["spectral checkpoint serialization", "restored clock/step/event schedule", "restored RNG probe (32 UInt64 draws)",
            "restored convection configuration", "every post-restart timestep/event schedule",
            "post-restart event locations and injection flags", "final clock/step/event schedule",
            "final RNG probe (32 UInt64 draws)"]
        open(joinpath(output, "report.toml.partial"), "w") do io
            TOML.print(io, report; sorted=true)
        end
        mv(joinpath(output, "report.toml.partial"), joinpath(output, "report.toml"))
        open(joinpath(output, "RESTART_VALIDATION_PASSED"), "w") do io
            println(io, "Same-device checkpoint equivalence passed on ", device)
        end
        println("RESTART VALIDATION PASSED on ", device, "; results: ", output)
        return report
    catch err
        open(joinpath(output, "FAILED.txt"), "w") do io
            showerror(io, err)
        end
        rethrow()
    end
end

function main(args)
    config, output, device = restart_options(args)
    validate_restart(config, output, device)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
