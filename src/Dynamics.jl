"One active layer; `prob.sol` contains s=(Laplacian-LD^-2)psi, not total PV."
struct Simulation{P}
    prob::P
    config::Dict{String,Any}
    x::Vector{Float64}
    y::Vector{Float64}
    background::Matrix{Float64}
    psid::Matrix{Float64}
    reference_s::Matrix{Float64}
    convection::ConvectionState
    device::String
    nu::Float64
    background_frequency_bound::Float64
end

smooth_window(r, ri, ro) = r <= ri ? 1.0 : r >= ro ? 0.0 :
    (t = (r-ri)/(ro-ri); 1-10t^3+15t^4-6t^5)

"Analytic integral of r*window(r); matches r²/2 in the science disk."
function planetary_pv(r, gamma, ri, ro)
    r <= ri && return -gamma*r^2/2
    d = ro-ri; t = min((r-ri)/d, 1.0)
    integral = ri^2/2 + ri*d*(t-2.5t^4+3t^5-t^6) + d^2*(t^2/2-2t^5+2.5t^6-6t^7/7)
    return -gamma*integral
end

# Smooth at the origin: the radial jet's streamfunction is analytic in r².
ring_psi(r, U, R, w) = U*w*(tanh((r^2-R^2)/(2R*w))-1)
cyclone_psi(r, U, R) = -U*R*exp(0.5-r^2/(2R^2))

function initial_streamfunctions(c, x, y)
    a, d, b = c["initial"], c["physics"], c["boundary"]
    matched = get(d, "deep_profile", "annular") == "matched_upper"
    psi = zeros(length(x), length(y)); psid = similar(psi); base = similar(psi)
    rng = MersenneTwister(c["experiment"]["seed"])
    modes = a["perturbation"] == "none" ? Int[] : a["perturbation"] == "mode" ? [a["mode"]] : collect(1:a["m_max"])
    phases = a["perturbation"] == "mode" ? [0.0] : 2pi .* rand(rng, length(modes))
    amp = isempty(modes) ? 0.0 : a["displacement_m"]/sqrt(length(modes))
    for j in eachindex(y), i in eachindex(x)
        r = hypot(x[i], y[j]); theta = atan(y[j], x[i])
        w = smooth_window(r, b["inner_radius_m"], b["outer_radius_m"])
        base[i,j] = w*(ring_psi(r, a["jet_speed_m_s"], a["jet_radius_m"], a["jet_width_m"]) + cyclone_psi(r, a["cyclone_speed_m_s"], a["cyclone_radius_m"]))
        psid[i,j] = matched ? 0.0 : w*ring_psi(r, d["deep_jet_speed_m_s"], d["deep_jet_radius_m"], d["deep_jet_width_m"])
        perturb = 0.0
        rho = r/a["jet_radius_m"]; z = (r^2-a["jet_radius_m"]^2)/(2a["jet_radius_m"]*a["jet_width_m"])
        for (m, phase) in zip(modes, phases)
            envelope = rho^m * exp(-m*(rho^2-1)/2)
            perturb -= amp*a["jet_speed_m_s"]*envelope/cosh(z)^2*cos(m*theta+phase)
        end
        psi[i,j] = base[i,j] + w*perturb
    end
    # Consistent zero-area-mean gauges are part of the convection definition.
    psi .-= mean(psi); base .-= mean(base)
    if matched
        # Prescribe the unperturbed ring AND cyclone, with the upper base's
        # taper and gauge. Never copy the initial nonaxisymmetric disturbance.
        psid .= get(d, "deep_match_fraction", 1.0) .* base
    else
        psid .-= mean(psid)
    end
    return psi, psid, base
end

function build_simulation(c::AbstractDict; device="cpu")
    validate_config(c)
    device in ("cpu", "gpu") || error("device must be cpu or gpu")
    if device == "gpu"
        CUDA.functional(true) || error("GPU requested but CUDA is not functional; CPU fallback is forbidden")
        CUDA.allowscalar(false)
    end
    dev = device == "gpu" ? GPU() : CPU()
    n = c["grid"]["nx"]; L = c["grid"]["length_m"]
    x = collect(range(-L/2, step=L/n, length=n)); y = copy(x)
    psi, psid, base = initial_streamfunctions(c, x, y)
    ld = c["physics"]["deformation_radius_m"]; kap2 = isinf(ld) ? 0.0 : 1/ld^2
    gamma = 2c["planet"]["rotation_rate_s"]/c["planet"]["radius_m"]^2
    b = c["boundary"]
    B = [kap2*psid[i,j] + planetary_pv(hypot(x[i], y[j]), gamma, b["inner_radius_m"], b["outer_radius_m"]) for i in 1:n, j in 1:n]
    B .-= mean(B)
    # Use the library's wave numbers and transform conventions everywhere.
    grid = TwoDGrid(CPU(); nx=n, ny=n, Lx=L, Ly=L, aliased_fraction=1/3,
                    nthreads=Threads.nthreads(), effort=FFTW.ESTIMATE, T=Float64)
    k2 = Array(grid.Krsq)
    project(a) = (ah = rfft(a); dealias!(ah, grid); irfft(ah, n))
    psid = project(psid)
    B = project(B)
    if get(c["physics"], "deep_profile", "annular") == "matched_upper"
        # Initial s is transformed/dealiased again by QG.set_q!. The pinned
        # library's cutoff mask is not exactly idempotent after a real-space
        # round trip at its conjugate cutoff edge. Apply the same two passes
        # to the prescribed base so matching remains exact on the stored grid.
        matched_psid = project(psid)
        B .+= kap2 .* (matched_psid .- psid)
        psid = matched_psid
    end
    s = project(irfft(-(k2 .+ kap2).*rfft(psi), n))
    sref = project(irfft(-(k2 .+ kap2).*rfft(base), n))
    taper = [1-smooth_window(hypot(xx,yy), b["inner_radius_m"], b["outer_radius_m"]) for xx in x, yy in y]
    kcut = 2pi/L*floor(n/3)
    p = c["numerics"]["hyper_order"]
    nu = 1/(c["numerics"]["cutoff_damping_time_s"]*kcut^(2p))
    # A deterministic optional boundary sponge; never a stochastic RK callback.
    da = device == "gpu" ? CUDA.CuArray : identity
    target = da(sref); mask = da(project(taper)); buffer = da(zeros(n,n))
    sigma = b["sponge_rate_s"]
    function sponge!(Fh, sol, t, clock, vars, params, g)
        if sigma == 0
            fill!(Fh, 0)
        else
            copyto!(vars.qh, sol)
            ldiv!(buffer, g.rfftplan, vars.qh)
            @. buffer = -sigma*mask*(buffer-target)
            buffer .-= sum(buffer)/length(buffer)
            mul!(Fh, g.rfftplan, buffer)
            dealias!(Fh, g)
        end
        return nothing
    end
    # QG.Problem does not expose FFTW's thread count. Assemble its pinned
    # constructors using our grid to honor the scheduler's Julia CPU allocation.
    modelgrid = device == "cpu" ? grid : TwoDGrid(dev; nx=n, ny=n, Lx=L, Ly=L,
        aliased_fraction=1/3, nthreads=Threads.nthreads(), effort=FFTW.ESTIMATE, T=Float64)
    params = isinf(ld) ? QG.BarotropicQGParams(modelgrid, 0.0, 0.0, B, 0.0, nu, p, sponge!) :
        QG.EquivalentBarotropicQGParams(modelgrid, ld, 0.0, 0.0, B, 0.0, nu, p, sponge!)
    vars = QG.ForcedVars(modelgrid)
    prob = QG.FourierFlows.Problem(QG.Equation(params, modelgrid), "RK4",
        c["numerics"]["dt_s"], modelgrid, vars, params)
    QG.set_q!(prob, da(s))
    conv = ConvectionState(c["convection"], c["experiment"]["seed"]+104729)
    Bh = rfft(B)
    gradB = maximum(hypot.(irfft(im .* grid.kr .* Bh, n), irfft(im .* grid.l .* Bh, n)))
    response = sqrt.(k2)./(k2 .+ kap2)
    response[1,1] = 0.0
    wavebound = gradB*maximum(response)
    return Simulation(prob, deepcopy(c), x, y, B, psid, sref, conv, device, nu, wavebound)
end

function fields(sim)
    QG.updatevars!(sim.prob)
    v = sim.prob.vars
    s, psi, u, vvel = Array(v.q), Array(v.ψ), Array(v.u), Array(v.v)
    return (; s, psi, u, v=vvel, totalpv=s.+sim.background)
end

function invariants(sim)
    f = fields(sim); ld = sim.config["physics"]["deformation_radius_m"]
    energy = mean(f.u.^2 .+ f.v.^2 .+ (f.psi./ld).^2)/2
    return (; energy, total_enstrophy=mean(abs2, f.totalpv)/2, mean_s=mean(f.s))
end

function check_step(sim, dt)
    p = sim.prob; QG.updatevars!(p)
    dx = sim.config["grid"]["length_m"]/length(sim.x)
    speed_x, speed_y = maximum(abs, p.vars.u), maximum(abs, p.vars.v)
    cfl = dt*(speed_x+speed_y)/dx
    isfinite(cfl) || error("nonfinite velocity at t=$(p.clock.t)")
    cfl <= sim.config["numerics"]["cfl_limit"] || error("CFL=$cfl exceeds limit; reduce dt_s (no silent adaptive change)")
    dt*sim.background_frequency_bound <= 0.5 || error("Background-PV wave bound exceeded; reduce dt_s")
    damping = dt*(sim.nu*maximum(p.grid.Krsq)^sim.config["numerics"]["hyper_order"] + sim.config["boundary"]["sponge_rate_s"])
    damping <= 2.0 || error("explicit RK4 damping bound exceeded ($damping); reduce dt_s or damping")
    return cfl
end

"Advance one step, landing exactly on the next convective event or `until`."
function advance!(sim; until=sim.config["numerics"]["duration_s"])
    p = sim.prob; t = p.clock.t
    until > t || return nothing
    dt = min(sim.config["numerics"]["dt_s"], until-t, sim.convection.next_event-t)
    if dt > 16eps(max(abs(t), 1.0))
        check_step(sim, dt)
        p.clock.dt = dt
        stepforward!(p, 1)
    else
        p.clock.t = min(until, sim.convection.next_event)
    end
    if p.clock.t >= sim.convection.next_event - 16eps(max(abs(p.clock.t), 1.0))
        p.clock.t = sim.convection.next_event
        f = fields(sim); before = invariants(sim)
        sampler = (xp,yp) -> local_pv_anomaly(f.totalpv, sim.x, sim.y, xp, yp)
        event = kick!(sim.convection, f.s, f.psi, sim.psid, f.totalpv, sim.x, sim.y, p.clock.t; qprime_sampler=sampler)
        QG.set_q!(p, sim.device == "gpu" ? CUDA.CuArray(f.s) : f.s)
        after = invariants(sim)
        return merge(event, (; delta_energy=after.energy-before.energy, delta_total_enstrophy=after.total_enstrophy-before.total_enstrophy))
    end
    return nothing
end
