using Test, SaturnPolarQG, Statistics, LinearAlgebra

function deep_profile_config(; alpha=1.0, ld=2e6, perturbation="none")
    c = load_config(joinpath(@__DIR__, "..", "configs", "smoke.toml"))
    c["convection"]["enabled"] = false
    c["initial"]["perturbation"] = perturbation
    c["physics"]["deformation_radius_m"] = ld
    c["physics"]["deep_profile"] = "matched_upper"
    c["physics"]["deep_match_fraction"] = alpha
    c["physics"]["deep_jet_speed_m_s"] = 0.0
    return c
end

@testset "Deep-profile contract and unchanged legacy default" begin
    legacy = load_config(joinpath(@__DIR__, "..", "configs", "smoke.toml"))
    original = deepcopy(legacy)
    @test validate_config(legacy) === legacy
    @test legacy == original
    explicit = deepcopy(legacy); explicit["physics"]["deep_profile"] = "annular"
    old, named = build_simulation(legacy), build_simulation(explicit)
    @test old.psid == named.psid
    @test old.background == named.background
    @test old.prob.sol == named.prob.sol
    d, b = legacy["physics"], legacy["boundary"]
    expected = [SaturnPolarQG.smooth_window(hypot(x,y),b["inner_radius_m"],b["outer_radius_m"])*
        SaturnPolarQG.ring_psi(hypot(x,y),d["deep_jet_speed_m_s"],d["deep_jet_radius_m"],d["deep_jet_width_m"])
        for x in old.x, y in old.y]
    expected .-= mean(expected)
    _, rawdeep, _ = SaturnPolarQG.initial_streamfunctions(legacy,old.x,old.y)
    @test rawdeep == expected
    for alpha in (0.0, 0.5, 1.0, 1.5)
        @test validate_config(deep_profile_config(; alpha)) isa AbstractDict
    end
    default = deep_profile_config(); delete!(default["physics"],"deep_match_fraction")
    @test build_simulation(default).psid == build_simulation(deep_profile_config()).psid
    for alpha in (-0.1, Inf, NaN, true, "1")
        bad = deep_profile_config(); bad["physics"]["deep_match_fraction"] = alpha
        @test_throws ErrorException validate_config(bad)
    end
    for keyvalue in (("deep_profile","other"),("deep_profile",true),("deep_jet_speed_m_s",5.0))
        bad = deep_profile_config(); bad["physics"][keyvalue[1]] = keyvalue[2]
        @test_throws ErrorException validate_config(bad)
    end
    bad = deepcopy(legacy); bad["physics"]["deep_match_fraction"] = 1.0
    @test_throws ErrorException validate_config(bad)
end

@testset "Full matching fixes upper velocity and cancels base stretching after projection" begin
    reference = build_simulation(deep_profile_config(; ld=Inf))
    fref = fields(reference)
    for ld in (0.5e6, 2e6, Inf), alpha in (0.0, 0.5, 1.0)
        sim = build_simulation(deep_profile_config(; ld, alpha))
        f = fields(sim)
        @test f.psi ≈ fref.psi rtol=5e-14
        @test f.u ≈ fref.u rtol=5e-14
        @test f.v ≈ fref.v rtol=5e-14
        @test sim.psid ≈ alpha .* reference.psid rtol=5e-14 atol=1e-9
        if alpha == 1
            @test f.psi ≈ sim.psid rtol=5e-14
            @test f.totalpv ≈ fref.totalpv rtol=5e-13
            p = SaturnPolarQG.physical_diagnostics(sim, f)
            if isfinite(ld)
                @test p.max_abs_hprime_over_H < 1e-13
                @test p.min_h_over_H ≈ 1.0 atol=1e-13
            else
                @test isnan(p.max_abs_hprime_over_H)
                @test isnan(p.min_h_over_H)
            end
            @test p.max_abs_zeta_over_f > 0
        end
    end
end

@testset "Matched deep flow includes cyclone and excludes initial disturbances" begin
    c = deep_profile_config()
    unperturbed = build_simulation(c)
    for disturbance in ("mode", "broadband")
        changed = deep_profile_config(; perturbation=disturbance)
        changed["experiment"]["seed"] = 127
        sim = build_simulation(changed)
        @test sim.psid == unperturbed.psid
        @test norm(fields(sim).psi-sim.psid) > 1e5
        @test sim.reference_s == unperturbed.reference_s
    end
    no_core = deepcopy(c); no_core["initial"]["cyclone_speed_m_s"] = 0.0
    @test norm(unperturbed.psid-build_simulation(no_core).psid) > 1e6
    psi, psid, base = SaturnPolarQG.initial_streamfunctions(c,unperturbed.x,unperturbed.y)
    @test psi == base == psid
    half = deepcopy(c); half["physics"]["deep_match_fraction"] = 0.5
    _, halfdeep, _ = SaturnPolarQG.initial_streamfunctions(half,unperturbed.x,unperturbed.y)
    @test halfdeep == 0.5 .* base
end

@testset "Thickness sign, local Coriolis, science disk, and CSV labels" begin
    sim = build_simulation(deep_profile_config())
    c = sim.config; f0 = 2c["planet"]["rotation_rate_s"]
    kap2 = inv(c["physics"]["deformation_radius_m"]^2)
    r2 = [x^2+y^2 for x in sim.x, y in sim.y]
    local_f = f0 .* (1 .- r2 ./ (2c["planet"]["radius_m"]^2))
    outside = r2 .> c["boundary"]["inner_radius_m"]^2
    psi = sim.psid .- 0.2*f0/kap2
    s = 0.3 .* local_f .- kap2 .* psi
    # Deliberately extreme exterior values must not pollute science-disk maxima.
    psi[outside] .= 1e16; s[outside] .= 1e6
    p = SaturnPolarQG.physical_diagnostics(sim, (; psi, s))
    @test p.max_abs_hprime_over_H ≈ 0.2 atol=1e-15
    @test p.min_h_over_H ≈ 0.8 atol=1e-15
    @test p.max_abs_zeta_over_f ≈ 0.3 atol=1e-15
    d = SaturnPolarQG.diagnostics(sim)
    @test d.header[9:11] == ["max_abs_hprime_over_H","min_h_over_H","max_abs_zeta_over_f"]
    @test d.header[12] == "m1_amplitude_m"
    @test length(d.header) == length(d.vals)
end

@testset "Matched profile stochastic restart and configuration identity" begin
    c = deep_profile_config(; perturbation="mode")
    c["convection"]["enabled"] = true
    whole, first = build_simulation(c), build_simulation(c)
    while whole.prob.clock.t < c["numerics"]["duration_s"]
        advance!(whole)
    end
    while first.prob.clock.t < 180.0
        advance!(first)
    end
    mktempdir() do dir
        path = joinpath(dir,"checkpoint.jld2")
        save_checkpoint(path, first)
        resumed = restore_checkpoint!(build_simulation(c),path)
        while resumed.prob.clock.t < c["numerics"]["duration_s"]
            advance!(resumed)
        end
        @test resumed.psid == whole.psid
        @test resumed.prob.clock.t == whole.prob.clock.t
        @test resumed.convection.events == whole.convection.events >= 4
        @test resumed.prob.sol ≈ whole.prob.sol rtol=1e-13
        bad = deepcopy(c); bad["physics"]["deep_match_fraction"] = 0.5
        @test_throws ErrorException restore_checkpoint!(build_simulation(bad),path)
    end
end
