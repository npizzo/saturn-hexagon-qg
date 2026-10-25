module Convection

using Random
using Statistics

export ConvectionState, kick!

"""
State for discrete, Thomson–McIntyre-inspired convective PV-pair injections.

The first event is scheduled from time zero. Save/restore this entire object,
including `rng`, to resume the same event sequence. No randomness belongs in
the time integrator's right-hand side. See `docs/convection.md` for the model
adaptations and configuration units.
"""
mutable struct ConvectionState
    next_event::Float64
    events::Int
    rng::MersenneTwister
    cfg::Dict{String,Any}
end

function ConvectionState(cfg::AbstractDict, seed::Integer)
    c = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in cfg)
    c["enabled"] = get(c, "enabled", false)
    c["enabled"] isa Bool || throw(ArgumentError("convection.enabled must be a Boolean"))
    rng = MersenneTwister(seed)
    if !c["enabled"]
        return ConvectionState(Inf, 0, rng, c)
    end
    for k in ("r_min_m", "r_max_m", "vortex_radius_m", "q_max_s", "bias_max",
              "psi_limit_m2_s", "interval_min_s", "interval_max_s")
        haskey(c, k) || throw(ArgumentError("missing convection.$k"))
        c[k] = Float64(c[k])
        isfinite(c[k]) || throw(ArgumentError("convection.$k must be finite"))
    end
    0 <= c["r_min_m"] < c["r_max_m"] || throw(ArgumentError("convection requires 0 <= r_min_m < r_max_m"))
    for k in ("vortex_radius_m", "q_max_s", "psi_limit_m2_s", "interval_min_s")
        c[k] > 0 || throw(ArgumentError("convection.$k must be positive"))
    end
    0 <= c["bias_max"] <= 1 || throw(ArgumentError("convection.bias_max must lie in [0, 1]"))
    c["interval_max_s"] >= c["interval_min_s"] || throw(ArgumentError("interval_max_s must be >= interval_min_s"))
    state = ConvectionState(0.0, 0, rng, c)
    state.next_event = _interval(state)
    return state
end

_interval(state) = state.cfg["interval_min_s"] + rand(state.rng) *
    (state.cfg["interval_max_s"] - state.cfg["interval_min_s"])

function _spacing(z::Vector{Float64}, label)
    length(z) >= 2 || throw(ArgumentError("$label needs at least two grid points"))
    all(isfinite, z) || throw(ArgumentError("$label must be finite"))
    d = z[2] - z[1]
    d > 0 || throw(ArgumentError("$label must be increasing"))
    all(v -> isapprox(v, d; rtol=1e-10, atol=1e-12 * abs(d)), diff(z)) ||
        throw(ArgumentError("$label must be uniformly spaced and omit its periodic endpoint"))
    return d
end

"Periodic bilinear interpolation, with array dimensions ordered as (x, y)."
function _sample(a, xp, yp, x, y, dx, dy)
    nx, ny = size(a)
    fx = mod((xp - x[1]) / dx, nx)
    fy = mod((yp - y[1]) / dy, ny)
    i0, j0 = floor(Int, fx), floor(Int, fy)
    ax, ay = fx - i0, fy - j0
    i, j = i0 + 1, j0 + 1
    ip, jp = mod1(i + 1, nx), mod1(j + 1, ny)
    return (1-ay) * ((1-ax) * a[i,j] + ax * a[ip,j]) +
           ay * ((1-ax) * a[i,jp] + ax * a[ip,jp])
end

_periodic_distance(d, length) = mod(d + length/2, length) - length/2

# The numerical profile is compact and parabolic, as in TM16 equation (3.1).
_profile(d2, r2) = max(0.0, 1.0 - d2 / r2)

function _strength(coldness, qprime, qmax, bmax)
    ramp = clamp(0.5 + coldness, 0.0, 1.0)
    bias = bmax * ramp
    ratio = 1.0 - bias
    # Once the background itself violates the cap there is no admissible
    # positive injection under this rule. Never allow a negative amplitude.
    if abs(qprime) >= qmax
        return (0.0, 0.0, bias)
    end
    a = min(qmax * ramp, qmax + qprime)
    if ratio > 0
        a = min(a, (qmax - qprime) / ratio)
    end
    a = max(0.0, a)
    return (ratio * a, -a, bias)
end

function _info(; attempted=false, injected=false, event=0, time_s=NaN,
               next_event_s=Inf, x_m=NaN, y_m=NaN, radius_m=NaN,
               cyclone_x_m=NaN, cyclone_y_m=NaN,
               anticyclone_x_m=NaN, anticyclone_y_m=NaN,
               q_cyclone_s=0.0, q_anticyclone_s=0.0,
               coldness=NaN, local_qprime_s=NaN, bias=NaN,
               removed_mean_s=0.0, increment_mean_s=0.0,
               increment_rms_s=0.0, increment_max_s=0.0,
               increment_min_s=0.0)
    return (; attempted, injected, event, time_s, next_event_s, x_m, y_m,
            radius_m, cyclone_x_m, cyclone_y_m, anticyclone_x_m,
            anticyclone_y_m, q_cyclone_s, q_anticyclone_s, coldness,
            local_qprime_s, bias, removed_mean_s, increment_mean_s,
            increment_rms_s, increment_max_s, increment_min_s)
end

"""
    kick!(state, s, psi, psid, qprime, x, y, time; qprime_sampler=nothing) -> NamedTuple

Apply one due impulsive vortex pair to `s`, sampling its midpoint uniformly by
area in the configured annulus. All arrays are CPU `Matrix{Float64}` with
dimensions `(length(x), length(y))`; coordinate vectors omit the periodic end
point. `qprime` is the instantaneous departure of total PV from its azimuthal
mean, computed by the caller. `psi` and `psid` must use consistent gauges.
If `qprime_sampler(x,y)` is supplied, it replaces interpolation of the qprime
matrix, allowing a direct circular-mean estimate at the randomly chosen point.

The caller must land on `state.next_event`; this routine intentionally applies
at most one event per call. Early/disabled calls do nothing. A scheduled event
with zero amplitude still increments the event counter and schedules the next
one. Returned strengths precede domain-mean compensation. Energy budgets must
be evaluated by the caller using its own inversion operator.
"""
function kick!(state::ConvectionState, s::Matrix{Float64},
               psi::Matrix{Float64}, psid::Matrix{Float64},
               qprime::Matrix{Float64}, x::Vector{Float64},
               y::Vector{Float64}, time::Real; qprime_sampler=nothing)
    t = Float64(time)
    isfinite(t) || throw(ArgumentError("event time must be finite"))
    if !state.cfg["enabled"] || t < state.next_event - 8eps(max(abs(t), 1.0))
        return _info(event=state.events, time_s=t, next_event_s=state.next_event)
    end
    shape = (length(x), length(y))
    all(a -> size(a) == shape, (s, psi, psid, qprime)) ||
        throw(DimensionMismatch("all convection arrays must have dimensions (length(x), length(y))"))
    dx, dy = _spacing(x, "x"), _spacing(y, "y")
    c = state.cfg
    r0 = c["vortex_radius_m"]
    # Ensure the two disks do not overlap. In the polar adaptation their
    # connecting line follows the local tangent at the sampled midpoint.
    halfsep = r0 + max(dx, dy)/2
    lx, ly = length(x) * dx, length(y) * dy
    2halfsep + 2r0 < min(lx, ly) ||
        throw(ArgumentError("vortex pair is too large for the periodic domain"))
    radius = sqrt(c["r_min_m"]^2 + rand(state.rng) *
                  (c["r_max_m"]^2 - c["r_min_m"]^2))
    theta = 2pi * rand(state.rng)
    xp, yp = radius*cos(theta), radius*sin(theta)
    tangent_x, tangent_y = -sin(theta), cos(theta)
    order = iseven(state.events) ? 1.0 : -1.0
    cx, cy = xp + order*halfsep*tangent_x, yp + order*halfsep*tangent_y
    ax, ay = xp - order*halfsep*tangent_x, yp - order*halfsep*tangent_y
    coldness = (_sample(psid, xp, yp, x, y, dx, dy) -
                _sample(psi, xp, yp, x, y, dx, dy)) / c["psi_limit_m2_s"]
    qp = qprime_sampler === nothing ? _sample(qprime, xp, yp, x, y, dx, dy) : qprime_sampler(xp, yp)
    isfinite(coldness) && isfinite(qp) || throw(ArgumentError("sampled convection fields must be finite"))
    qc, qa, bias = _strength(coldness, qp, c["q_max_s"], c["bias_max"])
    state.events += 1
    state.next_event = t + _interval(state)

    common = (; attempted=true, event=state.events, time_s=t,
               next_event_s=state.next_event, x_m=xp, y_m=yp, radius_m=radius,
               cyclone_x_m=cx, cyclone_y_m=cy, anticyclone_x_m=ax,
               anticyclone_y_m=ay, q_cyclone_s=qc, q_anticyclone_s=qa,
               coldness, local_qprime_s=qp, bias)
    if qc == 0 && qa == 0
        return _info(; common...)
    end

    delta = similar(s)
    r2 = r0^2
    @inbounds for j in eachindex(y), i in eachindex(x)
        dc2 = _periodic_distance(x[i]-cx, lx)^2 + _periodic_distance(y[j]-cy, ly)^2
        da2 = _periodic_distance(x[i]-ax, lx)^2 + _periodic_distance(y[j]-ay, ly)^2
        delta[i,j] = qc*_profile(dc2, r2) + qa*_profile(da2, r2)
    end
    removed = mean(delta)
    delta .-= removed
    # Remove the remaining floating-point sum error as well. This correction
    # is at roundoff level and does not enforce a zero mean on preexisting s.
    delta[end] -= sum(delta)
    s .+= delta
    return _info(; common..., injected=any(!iszero, delta),
                 removed_mean_s=removed, increment_mean_s=mean(delta),
                 increment_rms_s=sqrt(mean(abs2, delta)),
                 increment_max_s=maximum(delta), increment_min_s=minimum(delta))
end

end # module Convection
