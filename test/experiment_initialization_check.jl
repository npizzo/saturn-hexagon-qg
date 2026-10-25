using SaturnPolarQG, TOML, LinearAlgebra, Test, Printf

# Short CPU integration check of the actual six prepared 1024² inputs.
# This does not evolve the scientific experiments or validate their stability.
function check_initial_states(stage)
    reference_u=nothing; reference_v=nothing
    @testset "Prepared screen initialization and one step" begin
        for line in readlines(joinpath(stage,"members.tsv"))
            _,name=split(line,'\t')
            cfg=load_config(joinpath(stage,name))
            sim=build_simulation(cfg;device="cpu")
            f=fields(sim)
            if reference_u===nothing
                reference_u=copy(f.u); reference_v=copy(f.v)
            end
            velocity_error=hypot(norm(f.u-reference_u),norm(f.v-reference_v))/hypot(norm(reference_u),norm(reference_v))
            @test velocity_error < 1e-10
            d=SaturnPolarQG.diagnostics(sim)
            @test d.polar.coverage==1.0
            @test all(isfinite,d.polar.mode_amplitudes)
            cfl=SaturnPolarQG.check_step(sim,cfg["numerics"]["dt_s"])
            advance!(sim;until=cfg["numerics"]["dt_s"])
            after=fields(sim)
            @test all(all(isfinite,getfield(after,k)) for k in (:psi,:s,:u,:v))
            @test sim.prob.clock.t==cfg["numerics"]["dt_s"]
            @printf("%s: initial wind relative L2 %.3e; CFL %.4f; initial A6 %.3f m; t %.1f s\n",
                name,velocity_error,cfl,d.polar.mode_amplitudes[6],sim.prob.clock.t)
            flush(stdout)
            sim=nothing; f=nothing; after=nothing; GC.gc()
        end
    end
end
check_initial_states(isempty(ARGS) ? joinpath(@__DIR__,"..","configs","experiments_01","screen") : ARGS[1])
