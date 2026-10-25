using Test, SaturnPolarQG, LinearAlgebra, Statistics

const RSP = SaturnPolarQG

function rostami_test_config(; ld=3.2e6, jet_fraction=0.0, core_fraction=0.0, perturbation="none")
    c = load_config(joinpath(@__DIR__,"..","configs","smoke.toml"))
    c["grid"]["nx"] = 128; c["grid"]["length_m"] = 40e6
    a = c["initial"]
    a["profile"] = "rostami2017"
    a["jet_speed_m_s"] = 81.92; a["jet_radius_m"] = 10.784e6; a["jet_width_m"] = 3.2e6/sqrt(3)
    a["cyclone_speed_m_s"] = 153.6; a["cyclone_radius_m"] = 3.2e6*(0.42/1.3)^(1/1.3)
    a["core_regularization_m"] = 500e3
    a["perturbation"] = perturbation; a["displacement_m"] = 1e3
    a["m_min"] = 2; a["m_max"] = 8
    d = c["physics"]
    d["deformation_radius_m"] = ld; d["deep_profile"] = "matched_components"
    d["deep_jet_speed_m_s"] = 0.0
    d["deep_jet_fraction"] = jet_fraction; d["deep_cyclone_fraction"] = core_fraction
    c["boundary"]["inner_radius_m"] = 16e6; c["boundary"]["outer_radius_m"] = 19e6
    c["output"]["jet_search_halfwidth_m"] = 3e6
    c["convection"]["enabled"] = false
    return c
end

@testset "Rostami dimensional profiles and pole regularity" begin
    a = rostami_test_config()["initial"]
    p = RSP.rostami_profile(a)
    @test p.cyclone_scale ≈ 3.2e6 rtol=1e-15
    @test RSP.rostami_velocity(p.jet_peak_radius,p;component=:jet) ≈ a["jet_speed_m_s"] rtol=2e-15
    @test RSP.rostami_velocity(p.cyclone_peak_radius,p;component=:cyclone) ≈ a["cyclone_speed_m_s"] rtol=2e-15
    @test p.jet_peak_radius ≈ a["jet_radius_m"] rtol=1e-13
    for component in (:jet,:cyclone,:total)
        u,ur,urr = RSP.rostami_velocity_derivatives(0.0,p;component)
        @test u == 0.0
        @test isfinite(ur)
        @test urr == 0.0
        @test RSP.rostami_velocity(1.0,p;component) ≈ ur rtol=1e-9
        for r in (200e3,1.5e6,10.784e6,13e6)
            h = 100.0
            vm,v,vp = [RSP.rostami_velocity(s,p;component) for s in (r-h,r,r+h)]
            u,du,ddu = RSP.rostami_velocity_derivatives(r,p;component)
            @test u ≈ v
            @test (vp-vm)/(2h) ≈ du rtol=2e-6 atol=2e-12
            @test (vp-2v+vm)/h^2 ≈ ddu rtol=3e-5 atol=2e-15
        end
    end
    # The fit away from the regularized core approaches the paper's singular
    # fit, rather than changing the dimensional radius when LD changes.
    tiny = deepcopy(a); tiny["core_regularization_m"] = 1.0
    pt = RSP.rostami_profile(tiny)
    @test pt.cyclone_peak_radius ≈ a["cyclone_radius_m"] rtol=1e-11
    @test RSP.rostami_velocity(a["jet_radius_m"]+a["jet_width_m"],p;component=:jet) ≈ a["jet_speed_m_s"]/exp(1) rtol=1e-13
    @test_throws ArgumentError RSP.rostami_velocity(-1.0,p)
    @test_throws ArgumentError RSP.rostami_velocity(1.0,p;component=:other)
end

@testset "Radial integration recovers velocity and converges independently of grid" begin
    p = RSP.rostami_profile(rostami_test_config()["initial"])
    t = RSP.rostami_streamfunctions(p,19e6)
    refined = RSP.rostami_streamfunctions(p,19e6;points_per_scale=256)
    for component in (:jet,:cyclone,:total)
        @test RSP.radial_streamfunction(19e6,t;component) == 0.0
        @test RSP.radial_streamfunction(20e6,t;component) == 0.0
        for r in (0.0,12345.0,200e3,1.5e6,10.784e6,13e6,18.9e6)
            psi = RSP.radial_streamfunction(r,t;component)
            @test psi ≈ RSP.radial_streamfunction(r,refined;component) rtol=5e-11 atol=1e-4
            r == 0 && continue
            h = 100.0
            derivative = (RSP.radial_streamfunction(r+h,t;component)-RSP.radial_streamfunction(r-h,t;component))/(2h)
            @test derivative ≈ RSP.rostami_velocity(r,p;component) rtol=2e-6 atol=1e-6
        end
    end
end

@testset "Paper-profile configuration and resolved core guards" begin
    c = rostami_test_config()
    @test validate_config(c) === c
    for (section,key,value) in (("initial","profile","other"),("initial","core_regularization_m",0.0),
            ("initial","core_regularization_m",1.0),("initial","m_min",0),("initial","m_min",9),
            ("initial","m_min",true),("physics","deep_jet_fraction",-1.0),
            ("physics","deep_cyclone_fraction",Inf),("physics","deep_match_fraction",1.0),
            ("physics","deep_jet_speed_m_s",1.0))
        bad = deepcopy(c); bad[section][key] = value
        @test_throws ErrorException validate_config(bad)
    end
    for (section,key) in (("initial","core_regularization_m"),("physics","deep_jet_fraction"),("physics","deep_cyclone_fraction"))
        bad = deepcopy(c); delete!(bad[section],key)
        @test_throws ErrorException validate_config(bad)
    end
    bad = deepcopy(c); bad["initial"]["profile"] = "annular"
    @test_throws ErrorException validate_config(bad)
    bad = deepcopy(c); bad["experiment"]["purpose"] = "selection"
    @test_throws ErrorException validate_config(bad)
    resolved = deepcopy(bad); resolved["grid"]["nx"] = 512
    @test validate_config(resolved) === resolved
    for kind in ("annular","matched_upper")
        bad = deepcopy(c); bad["physics"]["deep_profile"] = kind
        @test_throws ErrorException validate_config(bad)
    end
end

@testset "Component isolation, LD independence, and projected full matching" begin
    c = rostami_test_config()
    sim = build_simulation(c)
    initial = fields(sim)
    for ld in (Inf,2e6,4e6), (aj,ac) in ((0.0,0.0),(1.0,0.0),(0.0,1.0),(1.0,1.0))
        s = build_simulation(rostami_test_config(;ld,jet_fraction=aj,core_fraction=ac))
        f = fields(s)
        @test f.psi ≈ initial.psi rtol=5e-14
        @test f.u ≈ initial.u rtol=5e-14
        @test f.v ≈ initial.v rtol=5e-14
        if aj == ac == 1.0
            @test f.psi ≈ s.psid rtol=5e-14
            if isfinite(ld)
                @test RSP.physical_diagnostics(s,f).max_abs_hprime_over_H < 1e-12
            end
        end
    end
    j = build_simulation(rostami_test_config(;jet_fraction=1.0))
    k = build_simulation(rostami_test_config(;core_fraction=1.0))
    both = build_simulation(rostami_test_config(;jet_fraction=1.0,core_fraction=1.0))
    @test j.psid+k.psid ≈ both.psid rtol=5e-14 atol=1e-7
    @test norm(j.psid-k.psid) > 1e6
    for perturbation in ("mode","broadband")
        cc = rostami_test_config(;jet_fraction=1.0,core_fraction=1.0,perturbation)
        pp = build_simulation(cc)
        @test pp.psid == both.psid
        @test norm(fields(pp).psi-pp.psid) > 1e4
        @test pp.reference_s == both.reference_s
    end
end

@testset "Broadband mode interval and displacement convention" begin
    c = rostami_test_config(;perturbation="broadband")
    # Direct evaluation on the jet circle isolates angular harmonics without
    # Cartesian interpolation or grid projection.
    R = c["initial"]["jet_radius_m"]
    n = 256
    theta = 2pi .* (0:n-1)./n
    psi,_,base = RSP.initial_streamfunctions(c,[0.0;R.*cos.(theta)],[0.0;R.*sin.(theta)])
    # Gauge removal changes the DC level; subtract the origin's gauge offset.
    signal = [psi[k+1,k+1]-base[k+1,k+1]-(psi[1,1]-base[1,1]) for k in 1:n]
    amplitudes = [2abs(sum(signal.*exp.(-im*m.*theta)))/n for m in 1:12]
    expected = c["initial"]["displacement_m"]*c["initial"]["jet_speed_m_s"]/sqrt(7)
    @test maximum(amplitudes[[1,9,10,11,12]]) < 1e-7*expected
    @test all(isapprox.(amplitudes[2:8],expected;rtol=1e-10))
    c["initial"]["m_min"] = 6; c["initial"]["m_max"] = 6
    a = RSP.initial_streamfunctions(c,[0.0,R],[0.0])
    @test all(isfinite,a[1])
    @test RSP.rostami_perturbation_envelope(0.0,R,c["initial"]["jet_width_m"],2) == 0.0
    @test RSP.rostami_perturbation_envelope(R,R,c["initial"]["jet_width_m"],6) == 1.0
end
