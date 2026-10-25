using SaturnPolarQG

function main(args)
    isempty(args) && error("Usage: julia --project=. scripts/run.jl CONFIG --output NEW_DIRECTORY [--device cpu|gpu] [--resume CHECKPOINT]")
    config = args[1]; opts = Dict{String,String}()
    iseven(length(args)-1) || error("Each option requires a value")
    for i in 2:2:length(args)
        args[i] in ("--output", "--device", "--resume") || error("Unknown option $(args[i])")
        haskey(opts, args[i]) && error("Duplicate option $(args[i])")
        opts[args[i]] = args[i+1]
    end
    haskey(opts, "--output") || error("--output is required")
    sim = run_experiment(config; output=opts["--output"], device=get(opts, "--device", "cpu"), resume=get(opts, "--resume", nothing))
    println("Completed ", sim.prob.clock.step, " steps to t=", sim.prob.clock.t, " s; ", sim.convection.events, " convection events. Results: ", abspath(opts["--output"]))
end
main(ARGS)
