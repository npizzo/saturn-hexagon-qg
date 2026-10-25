using Test
using Statistics
using LinearAlgebra
using FFTW
using SaturnPolarQG
import GeophysicalFlows

const DYNAMICS_QG = GeophysicalFlows.SingleLayerQG

function dynamics_problem(; ld=0.8, beta=0.0, dt=0.005, eta=nothing)
    # Dimensionless, periodic, deliberately resolved analytic test cases.
    return DYNAMICS_QG.Problem(GeophysicalFlows.CPU(); nx=64, ny=64,
        Lx=2pi, Ly=2pi, deformation_radius=ld, β=beta, eta=eta,
        U=0.0, μ=0.0, ν=0.0, dt=dt, stepper="RK4",
        aliased_fraction=1/3, T=Float64)
end

function dynamics_set_streamfunction!(prob, psi, ld)
    kappa2 = isinf(ld) ? 0.0 : inv(ld^2)
    s = irfft(-(Array(prob.grid.Krsq) .+ kappa2) .* rfft(psi), prob.grid.nx)
    DYNAMICS_QG.set_q!(prob, s)
    return nothing
end

function dynamics_invariants(prob, ld)
    DYNAMICS_QG.updatevars!(prob)
    q = Array(prob.vars.q)
    psi = Array(prob.vars.ψ)
    u, v = Array(prob.vars.u), Array(prob.vars.v)
    eta = Array(prob.params.eta)
    return (; energy=mean(u.^2 .+ v.^2 .+ (psi./ld).^2)/2,
            total_enstrophy=mean(abs2, q .+ eta)/2,
            anomaly_enstrophy=mean(abs2, q)/2,
            mean_s=mean(q))
end

@testset "Analytic Helmholtz inversion and velocity convention" begin
    for ld in (0.7, Inf)
        prob = dynamics_problem(; ld)
        x, y = Array(prob.grid.x), Array(prob.grid.y)
        amplitude, kx, ky = 0.02, 2, 3
        psi = [amplitude*cos(kx*xx+ky*yy) for xx in x, yy in y]
        kappa2 = isinf(ld) ? 0.0 : inv(ld^2)
        # The analytic anomaly is supplied directly, so this independently
        # checks the library's inversion rather than reusing its inverse.
        s = -(kx^2+ky^2+kappa2).*psi
        DYNAMICS_QG.set_q!(prob, s)
        DYNAMICS_QG.updatevars!(prob)
        expected_u = [amplitude*ky*sin(kx*xx+ky*yy) for xx in x, yy in y]
        expected_v = [-amplitude*kx*sin(kx*xx+ky*yy) for xx in x, yy in y]
        @test maximum(abs, Array(prob.vars.ψ) .- psi) < 2e-14
        @test maximum(abs, Array(prob.vars.u) .- expected_u) < 2e-14
        @test maximum(abs, Array(prob.vars.v) .- expected_v) < 2e-14
        @test maximum(abs, Array(prob.vars.q) .- s) < 2e-13
        expected_energy = amplitude^2*(kx^2+ky^2+kappa2)/4
        @test dynamics_invariants(prob,ld).energy ≈ expected_energy rtol=2e-13
    end
end

function dynamics_rossby_error(ld, dt; duration=2.0)
    beta, amplitude, kx, ky = 3.0, 0.02, 2, 1
    prob = dynamics_problem(; ld, beta, dt)
    x, y = Array(prob.grid.x), Array(prob.grid.y)
    psi0 = [amplitude*cos(kx*xx+ky*yy) for xx in x, yy in y]
    dynamics_set_streamfunction!(prob, psi0, ld)
    kappa2 = isinf(ld) ? 0.0 : inv(ld^2)
    omega = beta*kx/(kx^2+ky^2+kappa2)
    steps = round(Int,duration/dt)
    @assert isapprox(steps*dt,duration)
    GeophysicalFlows.stepforward!(prob,steps)
    DYNAMICS_QG.updatevars!(prob)
    actual = Array(prob.vars.ψ)
    # ∂t q + β∂x ψ = 0 gives cos(kx*x+ky*y+omega*t): westward phase.
    expected = [amplitude*cos(kx*xx+ky*yy+omega*duration) for xx in x, yy in y]
    phase = angle(fft(actual)[kx+1,ky+1]/fft(psi0)[kx+1,ky+1])
    return (; error=norm(actual.-expected)/norm(expected), phase,
            exact_phase=omega*duration, time=prob.clock.t)
end

@testset "Rossby-wave phase and fourth-order timestep convergence" begin
    for ld in (0.8, Inf)
        coarse = dynamics_rossby_error(ld,0.4)
        medium = dynamics_rossby_error(ld,0.2)
        fine = dynamics_rossby_error(ld,0.1)
        @test coarse.error > medium.error > fine.error > 0
        @test 12 < coarse.error/medium.error < 20
        @test 12 < medium.error/fine.error < 20
        @test fine.error < 5e-6
        @test fine.phase > 0
        @test fine.phase ≈ fine.exact_phase atol=5e-6
        @test fine.time ≈ 2.0 atol=5e-14
    end
end

@testset "Inviscid energy and total-PV enstrophy with fixed nonuniform PV" begin
    # Nonuniform eta makes anomaly-only enstrophy the wrong invariant. The
    # multimode flow also exercises nonlinear transfer, unlike a lone wave.
    ld = 0.8
    eta_fn(x,y) = 0.5cos(x+y) + 0.2sin(2x-y)
    prob = dynamics_problem(; ld, eta=eta_fn, dt=0.005)
    x, y = Array(prob.grid.x), Array(prob.grid.y)
    psi = [0.1cos(xx)+0.07sin(2yy)+0.03cos(2xx+yy) for xx in x, yy in y]
    dynamics_set_streamfunction!(prob,psi,ld)
    initial = dynamics_invariants(prob,ld)
    energy_errors, total_enstrophy_errors = Float64[], Float64[]
    last = initial
    for _ in 1:5
        GeophysicalFlows.stepforward!(prob,20)
        last = dynamics_invariants(prob,ld)
        push!(energy_errors,abs(last.energy/initial.energy-1))
        push!(total_enstrophy_errors,abs(last.total_enstrophy/initial.total_enstrophy-1))
        @test abs(last.mean_s-initial.mean_s) < 1e-14
    end
    @test maximum(energy_errors) < 1e-8
    @test maximum(total_enstrophy_errors) < 1e-9
    @test abs(last.anomaly_enstrophy/initial.anomaly_enstrophy-1) > 1e-6
end

@testset "Prescribed deep flow enters background PV, not initial upper velocity" begin
    smoke = load_config(joinpath(@__DIR__,"..","configs","smoke.toml"))
    cfg0 = deepcopy(smoke)
    cfg0["convection"]["enabled"] = false
    cfg0["physics"]["deformation_radius_m"] = 1.5e6
    cfg0["physics"]["deep_jet_speed_m_s"] = 0.0
    cfg1 = deepcopy(cfg0)
    cfg1["physics"]["deep_jet_speed_m_s"] = 30.0
    a, b = build_simulation(cfg0), build_simulation(cfg1)
    fa, fb = fields(a), fields(b)
    @test norm(b.psid-a.psid) > 0
    @test norm(b.background-a.background) > 0
    expected_delta = (b.psid-a.psid)/cfg0["physics"]["deformation_radius_m"]^2
    @test norm((b.background-a.background)-expected_delta) < 1e-12norm(expected_delta)
    @test fa.s ≈ fb.s rtol=2e-13
    @test fa.psi ≈ fb.psi rtol=2e-13
    @test fa.u ≈ fb.u rtol=2e-13
    @test fa.v ≈ fb.v rtol=2e-13
    @test norm(fb.totalpv-fa.totalpv) > 0

    cfg0["physics"]["deformation_radius_m"] = Inf
    cfg1["physics"]["deformation_radius_m"] = Inf
    ai, bi = build_simulation(cfg0), build_simulation(cfg1)
    @test ai.background == bi.background
    @test ai.prob.sol ≈ bi.prob.sol rtol=2e-13
    # With convection disabled the prescribed deep flow must disappear
    # completely from the barotropic-limit dynamics as well as its inversion.
    GeophysicalFlows.stepforward!(ai.prob,1)
    GeophysicalFlows.stepforward!(bi.prob,1)
    @test ai.prob.sol ≈ bi.prob.sol rtol=2e-13
    @test invariants(ai).energy ≈ invariants(bi).energy rtol=2e-13
end
