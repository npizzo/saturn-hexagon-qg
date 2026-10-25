using Test
using Statistics
using SaturnPolarQG.Convection

function convection_fixture(; kwargs...)
    cfg = Dict{String,Any}(
        "enabled"=>true, "r_min_m"=>8.0, "r_max_m"=>15.0,
        "vortex_radius_m"=>3.0, "q_max_s"=>2e-5, "bias_max"=>0.4,
        "psi_limit_m2_s"=>1.0, "interval_min_s"=>2.0, "interval_max_s"=>5.0)
    for (k,v) in kwargs
        cfg[String(k)] = v
    end
    x = collect(-32.0:1.0:31.0)
    y = copy(x)
    s = zeros(length(x), length(y))
    return cfg, s, zeros(size(s)), fill(1.0, size(s)), zeros(size(s)), x, y
end

@testset "Direct PV-anomaly sampling for the convection limiter" begin
    cfg, s, psi, psid, qp, x, y = convection_fixture()
    a, b = ConvectionState(cfg, 25), ConvectionState(cfg, 25)
    s2 = copy(s)
    e1 = kick!(a, s, psi, psid, qp, x, y, a.next_event)
    e2 = kick!(b, s2, psi, psid, fill(NaN,size(qp)), x, y, b.next_event; qprime_sampler=(xp,yp)->0.0)
    @test e1 == e2
    @test s == s2
    blocked = ConvectionState(cfg, 25)
    sb = zero(s)
    e3 = kick!(blocked, sb, psi, psid, qp, x, y, blocked.next_event; qprime_sampler=(xp,yp)->2cfg["q_max_s"])
    @test e3.attempted && !e3.injected
    @test iszero(sb)
end

@testset "Discrete convection configuration and event schedule" begin
    cfg, s, psi, psid, qp, x, y = convection_fixture()
    st = ConvectionState(cfg, 11)
    @test 2 <= st.next_event <= 5
    early = kick!(st, s, psi, psid, qp, x, y, 0.0)
    @test !early.attempted && st.events == 0 && iszero(s)
    t = st.next_event
    event = kick!(st, s, psi, psid, qp, x, y, t)
    @test event.attempted && event.injected && st.events == 1
    @test 2 <= st.next_event-t <= 5
    disabled = ConvectionState(Dict("enabled"=>false), 3)
    @test isinf(disabled.next_event)
    before = copy(s)
    @test !kick!(disabled, s, psi, psid, qp, x, y, 100.0).attempted
    @test s == before
    for (key, value) in (("r_min_m",-1), ("r_max_m",1), ("bias_max",1.1),
                         ("q_max_s",0), ("interval_max_s",1),
                         ("psi_limit_m2_s",NaN))
        invalid = copy(cfg); invalid[key] = value
        @test_throws ArgumentError ConvectionState(invalid, 1)
    end
end

@testset "Convection signs, bias, compensation, and tangent order" begin
    cfg, s, psi, psid, qp, x, y = convection_fixture()
    st = ConvectionState(cfg, 4)
    first = kick!(st, s, psi, psid, qp, x, y, st.next_event)
    @test first.q_cyclone_s > 0 > first.q_anticyclone_s
    @test first.q_cyclone_s / abs(first.q_anticyclone_s) ≈ 0.6
    @test maximum(s) > 0 > minimum(s)
    @test abs(sum(s)) < 1e-18
    @test first.removed_mean_s < 0
    @test first.increment_rms_s > 0
    tx, ty = -first.y_m/first.radius_m, first.x_m/first.radius_m
    @test (first.cyclone_x_m-first.anticyclone_x_m)*tx +
          (first.cyclone_y_m-first.anticyclone_y_m)*ty > 0
    second = kick!(st, s, psi, psid, qp, x, y, st.next_event)
    tx, ty = -second.y_m/second.radius_m, second.x_m/second.radius_m
    @test (second.cyclone_x_m-second.anticyclone_x_m)*tx +
          (second.cyclone_y_m-second.anticyclone_y_m)*ty < 0

    cfg["bias_max"] = 0.0
    unbiased = ConvectionState(cfg, 4)
    fill!(s, 7.0)  # compensation preserves any preexisting domain mean
    ev = kick!(unbiased, s, psi, psid, qp, x, y, unbiased.next_event)
    @test ev.q_cyclone_s == -ev.q_anticyclone_s
    @test mean(s) ≈ 7.0 atol=1e-14
    cfg["bias_max"] = 1.0
    anti = ConvectionState(cfg, 4)
    ev = kick!(anti, s, psi, psid, qp, x, y, anti.next_event)
    @test ev.q_cyclone_s == 0 && ev.q_anticyclone_s < 0
end

@testset "Coldness ramp and PV amplitude cap" begin
    cfg, s, psi, psid, qp, x, y = convection_fixture(bias_max=0.0)
    for (cold, fraction) in ((-1.0,0.0),(-0.5,0.0),(0.0,0.5),(0.5,1.0),(2.0,1.0))
        fill!(psid,cold); fill!(s,0)
        st = ConvectionState(cfg, 19)
        ev = kick!(st,s,psi,psid,qp,x,y,st.next_event)
        @test ev.q_cyclone_s ≈ cfg["q_max_s"]*fraction
        @test ev.q_anticyclone_s ≈ -cfg["q_max_s"]*fraction
        @test st.events == 1
    end
    fill!(psid,1.0)
    for background in (-4e-5,-2e-5,-1.9e-5,1.9e-5,2e-5,4e-5)
        fill!(qp,background); fill!(s,0)
        st = ConvectionState(cfg,19)
        ev = kick!(st,s,psi,psid,qp,x,y,st.next_event)
        @test ev.q_cyclone_s >= 0 && ev.q_anticyclone_s <= 0
        if abs(background) >= cfg["q_max_s"]
            @test !ev.injected && iszero(s)
        else
            @test max(ev.q_cyclone_s+background,-ev.q_anticyclone_s-background) <= cfg["q_max_s"]+1e-19
        end
    end
end

@testset "Seed and serialized-state-compatible restart reproducibility" begin
    cfg, s, psi, psid, qp, x, y = convection_fixture()
    a, b = ConvectionState(cfg,99), ConvectionState(cfg,99)
    sa, sb = copy(s), copy(s)
    for _ in 1:5
        ea = kick!(a,sa,psi,psid,qp,x,y,a.next_event)
        eb = kick!(b,sb,psi,psid,qp,x,y,b.next_event)
        @test ea == eb && sa == sb
    end
    resumed = deepcopy(a)
    sr = copy(sa)
    for _ in 1:8
        @test kick!(a,sa,psi,psid,qp,x,y,a.next_event) ==
              kick!(resumed,sr,psi,psid,qp,x,y,resumed.next_event)
        @test sa == sr
    end
end

@testset "Area-uniform midpoint sampling and periodic interpolation" begin
    cfg, s, psi, psid, qp, x, y = convection_fixture()
    fill!(psid,-1.0) # suppress amplitude; the event location still samples normally
    st = ConvectionState(cfg,121)
    radii, cosines, sines = Float64[], Float64[], Float64[]
    for _ in 1:1000
        ev = kick!(st,s,psi,psid,qp,x,y,st.next_event)
        @test cfg["r_min_m"] <= ev.radius_m <= cfg["r_max_m"]
        @test hypot(ev.x_m,ev.y_m) ≈ ev.radius_m
        push!(radii,(ev.radius_m^2-cfg["r_min_m"]^2)/(cfg["r_max_m"]^2-cfg["r_min_m"]^2))
        push!(cosines,ev.x_m/ev.radius_m); push!(sines,ev.y_m/ev.radius_m)
    end
    @test abs(mean(radii)-0.5) < 0.04
    @test abs(mean(cosines)) < 0.08 && abs(mean(sines)) < 0.08
    a = reshape(collect(1.0:16.0),4,4)
    xx = collect(0.0:1.0:3.0)
    sample = Convection._sample
    @test sample(a,0.5,0.5,xx,xx,1.0,1.0) == (a[1,1]+a[2,1]+a[1,2]+a[2,2])/4
    @test sample(a,3.5,0.5,xx,xx,1.0,1.0) == (a[4,1]+a[1,1]+a[4,2]+a[1,2])/4
    @test sample(a,7.5,-3.5,xx,xx,1.0,1.0) == sample(a,3.5,0.5,xx,xx,1.0,1.0)
end
