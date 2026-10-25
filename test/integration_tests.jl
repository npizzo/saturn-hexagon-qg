using Test, SaturnPolarQG, TOML, JLD2, LinearAlgebra

@testset "Integrated stochastic restart and volume constraint" begin
    c = load_config(joinpath(@__DIR__, "..", "configs", "smoke.toml"))
    whole = build_simulation(c)
    while whole.prob.clock.t < c["numerics"]["duration_s"]
        advance!(whole)
    end
    firsthalf = build_simulation(c)
    while firsthalf.prob.clock.t < 180.0
        advance!(firsthalf)
    end
    mktempdir() do dir
        checkpoint = joinpath(dir, "checkpoint.jld2")
        save_checkpoint(checkpoint, firsthalf)
        resumed = restore_checkpoint!(build_simulation(c), checkpoint)
        while resumed.prob.clock.t < c["numerics"]["duration_s"]
            advance!(resumed)
        end
        @test resumed.prob.clock.t == whole.prob.clock.t
        @test resumed.prob.clock.step == whole.prob.clock.step
        @test resumed.convection.events == whole.convection.events >= 4
        @test resumed.convection.next_event == whole.convection.next_event
        @test rand(deepcopy(resumed.convection.rng)) == rand(deepcopy(whole.convection.rng))
        @test norm(Array(resumed.prob.sol)-Array(whole.prob.sol))/norm(Array(whole.prob.sol)) < 1e-13
        @test abs(invariants(resumed).mean_s) < 1e-18
        changed = deepcopy(c); changed["physics"]["deformation_radius_m"] *= 1.1
        @test_throws ErrorException restore_checkpoint!(build_simulation(changed), checkpoint)
    end
end

@testset "Run artifacts and no-overwrite behavior" begin
    path = joinpath(@__DIR__, "..", "configs", "smoke.toml")
    mktempdir() do dir
        out = joinpath(dir, "run")
        sim = run_experiment(path; output=out)
        @test isfile(joinpath(out, "COMPLETED"))
        @test isfile(joinpath(out, "events.csv"))
        @test isfile(joinpath(out, "domain.jld2"))
        @test isfile(joinpath(out, "checkpoint.jld2"))
        @test TOML.parsefile(joinpath(out, "metadata.toml"))["device"] == "cpu"
        @test length(readlines(joinpath(out, "diagnostics.csv"))) >= 5
        @test length(readlines(joinpath(out, "profile_latest.csv"))) == sim.config["output"]["nr"]+1
        imported = SaturnPolarQG.PolarDiagnostics.load_wind_profile(joinpath(out, "profile_latest.csv"))
        @test length(imported.radius_m) == sim.config["output"]["nr"]
        @test_throws ErrorException run_experiment(path; output=out)
    end
end

@testset "Configuration and boundary science guards" begin
    c = load_config(joinpath(@__DIR__, "..", "configs", "smoke.toml"))
    bad = deepcopy(c); bad["convection"]["vortex_radius_m"] = 1e5
    @test_throws ErrorException validate_config(bad)
    bad = deepcopy(c); bad["numerics"]["dt_s"] = -1
    @test_throws ErrorException validate_config(bad)
    bad = deepcopy(c); bad["initial"]["mystery_parameter"] = 1
    @test_throws ErrorException validate_config(bad)
    gamma, ri, ro = 1e-19, 2e7, 2.5e7
    @test SaturnPolarQG.planetary_pv(1e7, gamma, ri, ro) == -gamma*1e14/2
    @test SaturnPolarQG.planetary_pv(3e7, gamma, ri, ro) == SaturnPolarQG.planetary_pv(4e7, gamma, ri, ro)
    h = 100.0
    slope = (SaturnPolarQG.planetary_pv(ri+h,gamma,ri,ro)-SaturnPolarQG.planetary_pv(ri-h,gamma,ri,ro))/(2h)
    @test slope ≈ -gamma*ri rtol=1e-8
    fast = deepcopy(c); fast["numerics"]["dt_s"] = 1e7
    @test_throws ErrorException SaturnPolarQG.check_step(build_simulation(fast), 1e7)
end
