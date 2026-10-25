"Dimensional, pole-regularized counterpart of the velocity fit in Rostami et al. (2017)."
struct RostamiProfile
    jet_speed::Float64
    jet_radius::Float64
    jet_width::Float64
    cyclone_speed::Float64
    cyclone_scale::Float64
    core_radius::Float64
    jet_peak_radius::Float64
    cyclone_peak_radius::Float64
    jet_normalization::Float64
    cyclone_normalization::Float64
end

const ROSTAMI_ALPHA = 0.42
const ROSTAMI_BETA = 1.3

function positive_peak(derivative, scale)
    lo, hi = 0.0, Float64(scale)
    for _ in 1:100
        derivative(hi) <= 0 && break
        hi *= 2
    end
    derivative(hi) <= 0 || error("Could not bracket radial-profile maximum")
    for _ in 1:80
        mid = (lo+hi)/2
        derivative(mid) > 0 ? (lo=mid) : (hi=mid)
    end
    return (lo+hi)/2
end

function rostami_jet_shape(r, R, w)
    # Odd continuation gives an analytic Cartesian streamfunction at the pole.
    # expm1 avoids cancellation in the difference of the two Gaussians.
    em = exp(-((r-R)/w)^2)
    ep = exp(-((r+R)/w)^2)
    u = em*(-expm1(-4r*R/w^2))
    ur = (-2(r-R)*em+2(r+R)*ep)/w^2
    urr = ((4(r-R)^2/w^2-2)*em-(4(r+R)^2/w^2-2)*ep)/w^2
    return (u, ur, urr)
end

function rostami_core_shape(r, L, epsilon)
    a, b = ROSTAMI_ALPHA, ROSTAMI_BETA
    s = hypot(r, epsilon)
    A = s^(a-1)/L^a * exp(-(s/L)^b)
    g = (a-1)/s^2-b*s^(b-2)/L^b
    h = -2(a-1)/s^4-b*(b-2)*s^(b-4)/L^b
    return (r*A, A*(1+r^2*g), A*r*(3g+r^2*(g^2+h)))
end

function rostami_profile(initial::AbstractDict)
    R, w = Float64(initial["jet_radius_m"]), Float64(initial["jet_width_m"])
    epsilon = Float64(initial["core_regularization_m"])
    L = Float64(initial["cyclone_radius_m"])/(ROSTAMI_ALPHA/ROSTAMI_BETA)^(1/ROSTAMI_BETA)
    all(x -> isfinite(x) && x > 0, (R,w,epsilon,L)) || error("Rostami profile lengths must be finite and positive")
    # Remove positive exponential factors before finding maxima, so distant
    # underflowing tails cannot masquerade as derivative zeros.
    jp = positive_peak(r -> R-r*tanh(2r*R/w^2), max(R,w))
    core_peak_condition(r) = (s=hypot(r,epsilon);
        1+(ROSTAMI_ALPHA-1)*r^2/s^2-ROSTAMI_BETA*r^2*s^(ROSTAMI_BETA-2)/L^ROSTAMI_BETA)
    cp = positive_peak(core_peak_condition, max(L,epsilon))
    jet_norm, core_norm = inv(rostami_jet_shape(jp,R,w)[1]), inv(rostami_core_shape(cp,L,epsilon)[1])
    isfinite(jet_norm) && isfinite(core_norm) || error("Profile normalization is not representable; check profile length scales")
    return RostamiProfile(Float64(initial["jet_speed_m_s"]), R, w,
        Float64(initial["cyclone_speed_m_s"]), L, epsilon, jp, cp,
        jet_norm, core_norm)
end

"Return (U, dU/dr, d²U/dr²); component is :jet, :cyclone, or :total."
function rostami_velocity_derivatives(r::Real, p::RostamiProfile; component=:total)
    r >= 0 || throw(ArgumentError("radial coordinate must be nonnegative"))
    component in (:jet,:cyclone,:total) || throw(ArgumentError("unknown velocity component"))
    j = component == :cyclone ? (0.0,0.0,0.0) :
        (p.jet_speed*p.jet_normalization) .* rostami_jet_shape(r,p.jet_radius,p.jet_width)
    c = component == :jet ? (0.0,0.0,0.0) :
        (p.cyclone_speed*p.cyclone_normalization) .* rostami_core_shape(r,p.cyclone_scale,p.core_radius)
    return j .+ c
end

rostami_velocity(r, p::RostamiProfile; component=:total) =
    rostami_velocity_derivatives(r,p;component)[1]

"Streamfunction table, integrating each velocity component to zero at rmax."
struct RostamiRadialTable
    rmax::Float64
    dr::Float64
    jet::Matrix{Float64}
    cyclone::Matrix{Float64}
end

function streamfunction_x_derivatives(r, p, component)
    u, ur, _ = rostami_velocity_derivatives(r,p;component)
    if r > 0
        return (u/(2r), (r*ur-u)/(4r^3))
    end
    # x=r² makes each interpolation polynomial exactly smooth at the origin.
    # ψ_x(0)=U'(0)/2 and ψ_xx(0)=U'''(0)/12.
    if component == :jet
        R, w = p.jet_radius,p.jet_width
        second = p.jet_speed*p.jet_normalization*exp(-(R/w)^2)*(4R^3/(3w^6)-2R/w^4)
    else
        a,b,e,L = ROSTAMI_ALPHA,ROSTAMI_BETA,p.core_radius,p.cyclone_scale
        A = e^(a-1)/L^a*exp(-(e/L)^b)
        g = (a-1)/e^2-b*e^(b-2)/L^b
        second = p.cyclone_speed*p.cyclone_normalization*A*g/4
    end
    return (ur/2,second)
end

function radial_component_table(p, rmax, n, component)
    h = rmax/n
    r = collect(range(0.0,rmax;length=n+1))
    u = [rostami_velocity(ri,p;component) for ri in r]
    psi = zeros(n+1)
    # Composite Simpson integration, accumulated from the zero exterior gauge.
    for i in n:-1:1
        psi[i] = psi[i+1]-h*(u[i]+4rostami_velocity((r[i]+r[i+1])/2,p;component)+u[i+1])/6
    end
    derivatives = [streamfunction_x_derivatives(ri,p,component) for ri in r]
    coefficients = zeros(6,n)
    for i in 1:n
        dx = (2i-1)*h^2
        a0, a1, a2 = psi[i], dx*derivatives[i][1], dx^2*derivatives[i][2]/2
        D = psi[i+1]-a0-a1-a2
        E = dx*derivatives[i+1][1]-a1-2a2
        F = dx^2*derivatives[i+1][2]-2a2
        coefficients[:,i] .= (a0,a1,a2,10*D-4*E+F/2,-15*D+7*E-F,6*D-3*E+F/2)
    end
    return coefficients
end

function rostami_streamfunctions(p::RostamiProfile, rmax::Real; points_per_scale=128)
    isfinite(rmax) && rmax > 0 || throw(ArgumentError("table outer radius must be finite and positive"))
    points_per_scale isa Integer && points_per_scale >= 8 || throw(ArgumentError("points_per_scale must be an integer >=8"))
    n = max(8,ceil(Int,rmax/min(p.jet_width,p.core_radius,p.cyclone_scale)*points_per_scale))
    j = radial_component_table(p,Float64(rmax),n,:jet)
    c = radial_component_table(p,Float64(rmax),n,:cyclone)
    return RostamiRadialTable(Float64(rmax),Float64(rmax)/n,j,c)
end

function radial_streamfunction(r::Real, table::RostamiRadialTable; component=:total)
    r >= 0 || throw(ArgumentError("radial coordinate must be nonnegative"))
    component in (:jet,:cyclone,:total) || throw(ArgumentError("unknown streamfunction component"))
    r >= table.rmax && return 0.0
    i = min(floor(Int,r/table.dr)+1,size(table.jet,2))
    t = (r^2-((i-1)*table.dr)^2)/((2i-1)*table.dr^2)
    value(a) = evalpoly(t,(a[1,i],a[2,i],a[3,i],a[4,i],a[5,i],a[6,i]))
    return component == :jet ? value(table.jet) : component == :cyclone ? value(table.cyclone) : value(table.jet)+value(table.cyclone)
end

"Smooth Gaussian annular seed; displacement_m is its streamfunction scale/Ujet."
function rostami_perturbation_envelope(r, R, w, m)
    # Matching the Gaussian's curvature at R; dependence on r² removes the
    # tiny nonsmooth polar tail of exp(-(r-R)²/w²) times an angular harmonic.
    rho = r/R
    return rho^m*exp(-m*(rho^2-1)/2-((r^2-R^2)/(2R*w))^2)
end
