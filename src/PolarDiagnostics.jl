module PolarDiagnostics

using Statistics

export radial_anomaly, local_pv_anomaly, polar_diagnostics, load_wind_profile

function _grid(x, y, dims)
    length(x) >= 2 && length(y) >= 2 || throw(ArgumentError("At least two grid points per direction are required"))
    dims == (length(x), length(y)) || throw(DimensionMismatch("Fields must have size (length(x), length(y))"))
    all(isfinite, x) && all(isfinite, y) || throw(ArgumentError("Coordinates must be finite"))
    dx, dy = Float64(x[2] - x[1]), Float64(y[2] - y[1])
    dx > 0 && dy > 0 || throw(ArgumentError("Coordinates must increase strictly"))
    all(d -> isapprox(d, dx; rtol=1e-8, atol=1e-10 * dx), diff(x)) || throw(ArgumentError("x must be uniformly spaced"))
    all(d -> isapprox(d, dy; rtol=1e-8, atol=1e-10 * dy), diff(y)) || throw(ArgumentError("y must be uniformly spaced"))
    return (; nx=length(x), ny=length(y), x0=Float64(first(x)), y0=Float64(first(y)), dx, dy)
end

# Coordinates describe a periodic, endpoint-excluded Cartesian grid. Masked
# corners are not filled: a sample touching a missing corner remains missing.
@inline function _sample(field, x, y, grid)
    fx = mod((x - grid.x0) / grid.dx, grid.nx)
    fy = mod((y - grid.y0) / grid.dy, grid.ny)
    ix, iy = floor(Int, fx), floor(Int, fy)
    ax, ay = fx - ix, fy - iy
    i, j = ix + 1, iy + 1
    ip, jp = mod1(i + 1, grid.nx), mod1(j + 1, grid.ny)
    f00, f10 = field[i, j], field[ip, j]
    f01, f11 = field[i, jp], field[ip, jp]
    all(isfinite, (f00, f10, f01, f11)) || return NaN
    return (1-ay) * ((1-ax)*f00 + ax*f10) + ay * ((1-ax)*f01 + ax*f11)
end

# Four-point Lagrange interpolation on nodes -1, 0, 1, 2 about the
# lower grid point. Tensor interpolation has fourth-order spatial accuracy
# for a smooth field and does not introduce the bias of radial bins.
@inline _cubic_weights(a) = (-a*(a-1)*(a-2)/6,
                             (a+1)*(a-1)*(a-2)/2,
                            -(a+1)*a*(a-2)/2,
                             (a+1)*a*(a-1)/6)

@inline function _sample_cubic(field, x, y, grid)
    fx = mod((x - grid.x0) / grid.dx, grid.nx)
    fy = mod((y - grid.y0) / grid.dy, grid.ny)
    ix, iy = floor(Int, fx), floor(Int, fy)
    wx, wy = _cubic_weights(fx-ix), _cubic_weights(fy-iy)
    base = Float64(field[ix+1,iy+1])
    isfinite(base) || return NaN
    delta = 0.0
    for b in 1:4, a in 1:4
        value = field[mod1(ix+a-1,grid.nx),mod1(iy+b-1,grid.ny)]
        isfinite(value) || return NaN
        # Interpolate differences to preserve a constant background exactly.
        delta += wx[a]*wy[b]*(value-base)
    end
    return base + delta
end

"""
    local_pv_anomaly(field, x, y, xp, yp; ntheta=nothing)

Evaluate `field(xp,yp) - mean_theta(field(r*cos(theta),r*sin(theta)))`,
where `r=hypot(xp,yp)`, using the same periodic tensor cubic Lagrange
interpolant for both terms. Arrays use `[x_index,y_index]`; the origin is
the pole. Default angular quadrature uses
`max(64, ceil(Int, 4π*r/min(dx,dy)))` equally spaced points, with no duplicated
endpoint. An explicit `ntheta` must be an integer of at least 8.

This pointwise definition is used by the convection limiter: it avoids the
radial-bin bias of `radial_anomaly`. Smooth axisymmetric fields still have a
small fourth-order Cartesian interpolation residual, so grid convergence
must be checked against the limiter's PV scale. The periodic sampler wraps
coordinates; physical applications should keep the whole circle inside
the usable polar cap. Any nonfinite value in a required 4×4 stencil yields
NaN; no missing sectors are filled. At the pole the anomaly is zero if the
sample is finite.
"""
function local_pv_anomaly(field::AbstractMatrix{<:Real}, x::AbstractVector,
                          y::AbstractVector, xp::Real, yp::Real; ntheta=nothing)
    grid = _grid(x,y,size(field))
    grid.nx >= 4 && grid.ny >= 4 || throw(ArgumentError("Cubic interpolation requires at least four points per direction"))
    isfinite(xp) && isfinite(yp) || throw(ArgumentError("Sample coordinates must be finite"))
    radius = hypot(xp,yp)
    n = if isnothing(ntheta)
        max(64,ceil(Int,4π*radius/min(grid.dx,grid.dy)))
    else
        ntheta isa Integer && ntheta >= 8 || throw(ArgumentError("ntheta must be an integer of at least 8"))
        Int(ntheta)
    end
    local_value = _sample_cubic(field,xp,yp,grid)
    isfinite(local_value) || return NaN
    radius == 0 && return 0.0
    # Accumulate differences, with compensation, to avoid subtracting two
    # large nearly identical PV backgrounds after taking the circle mean.
    total, compensation = 0.0, 0.0
    for k in 0:n-1
        st, ct = sincos(2π*k/n)
        value = _sample_cubic(field,radius*ct,radius*st,grid)
        isfinite(value) || return NaN
        term = (local_value-value) - compensation
        updated = total + term
        compensation = (updated-total)-term
        total = updated
    end
    return total/n
end

@inline function _linear_clamped(r, positions, values)
    r <= first(positions) && return first(values)
    r >= last(positions) && return last(values)
    i = searchsortedlast(positions, r)
    a = (r - positions[i]) / (positions[i+1] - positions[i])
    return values[i] + a*(values[i+1]-values[i])
end

"""
    radial_anomaly(field, x, y; nbins=max(16, min(length(x),length(y))÷2))

Subtract a radial-bin estimate of the azimuthal mean about (0,0). Array order
is `[x_index,y_index]`. Finite grid values in each bin are averaged; the bin's
mean radius supplies the interpolation coordinate. Empty bins are skipped and
the two end values are extended constantly. NaNs remain masked. This is an
approximate radial projection, not an exact spectral projector. It is retained
for grid diagnostics and is not used by the convection limiter; that uses
`local_pv_anomaly` to avoid radial-bin bias. Check binning and grid sensitivity
before interpreting this approximate array as a physical anomaly.
"""
function radial_anomaly(field::AbstractMatrix{<:Real}, x::AbstractVector, y::AbstractVector;
                        nbins::Integer=max(16, min(length(x), length(y)) ÷ 2))
    _grid(x, y, size(field))
    nbins >= 2 || throw(ArgumentError("nbins must be at least 2"))
    rmax = hypot(maximum(abs, x), maximum(abs, y))
    rmax > 0 || throw(ArgumentError("Grid must span nonzero radius"))
    sumr, sumf, counts = zeros(nbins), zeros(nbins), zeros(Int, nbins)
    for j in eachindex(y), i in eachindex(x)
        value = field[i,j]
        isfinite(value) || continue
        r = hypot(x[i], y[j])
        k = min(nbins, floor(Int, r / rmax * nbins) + 1)
        sumr[k] += r
        sumf[k] += value
        counts[k] += 1
    end
    keep = findall(>(0), counts)
    out = fill(NaN, size(field))
    isempty(keep) && return out
    positions, values = sumr[keep] ./ counts[keep], sumf[keep] ./ counts[keep]
    for j in eachindex(y), i in eachindex(x)
        isfinite(field[i,j]) || continue
        out[i,j] = field[i,j] - _linear_clamped(hypot(x[i], y[j]), positions, values)
    end
    return out
end

_finite_mean(a) = isempty(a) ? NaN : mean(a)

"""
    polar_diagnostics(u,v,totalpv,psi,s,x,y;
        jet_radius_m, jet_halfwidth_m, m_max=12, ntheta=360, nr=64,
        min_peak_contrast_fraction=1e-3)

Sample a north-polar annulus with a periodic bilinear interpolant. `u,v` are
Cartesian velocities in m/s with `u=-ψ_y, v=ψ_x`; positive azimuthal velocity
is counterclockwise. `s=(Δ-L_D^-2)ψ`; `totalpv=s+η`. The last two scalar
arrays are shape-checked for the common snapshot contract, but no energy is
computed here. CPU arrays have order `[x_index,y_index]`.

`jet_halfwidth_m` is the *search* halfwidth, not a fitted physical jet width.
The returned `radii` span this annulus. At each angle locate an interior,
positive azimuthal-speed maximum and refine it with a three-point parabola.
Flat, boundary-peaked, or missing profiles are invalid. `coverage` is the
valid angular fraction; `radial_coverage` separately reports each ring's
velocity coverage. Incomplete angular coverage produces NaN Fourier modes
rather than treating missing sectors as zero.
The radial speed range must exceed `min_peak_contrast_fraction` times its
maximum absolute speed (default 0.1%), to reject interpolation ripples in a
flat profile; this is a detection threshold, not an uncertainty estimate.

For displacement `A*cos(m*θ+φ)`, `mode_coefficients[m]=A/2*exp(im*φ)` and
`mode_amplitudes[m]=A`. `meanjet` is mean detected radius in meters;
`mean_peak_speed` is mean fitted peak speed in m/s. `eddy_uv` is the
azimuthal mean of `(u_r-mean(u_r))*(u_θ-mean(u_θ))` in m²/s².
"""
function polar_diagnostics(u::AbstractMatrix, v::AbstractMatrix, totalpv::AbstractMatrix,
                           psi::AbstractMatrix, s::AbstractMatrix,
                           x::AbstractVector, y::AbstractVector;
                           jet_radius_m::Real, jet_halfwidth_m::Real,
                           m_max::Integer=12, ntheta::Integer=360, nr::Integer=64,
                           min_peak_contrast_fraction::Real=1e-3)
    grid = _grid(x, y, size(u))
    all(a -> size(a) == size(u), (v, totalpv, psi, s)) || throw(DimensionMismatch("Snapshot fields must have identical shapes"))
    isfinite(jet_radius_m) && jet_radius_m > 0 || throw(ArgumentError("jet_radius_m must be finite and positive"))
    isfinite(jet_halfwidth_m) && jet_halfwidth_m > 0 || throw(ArgumentError("jet_halfwidth_m must be finite and positive"))
    ntheta >= 8 && nr >= 5 || throw(ArgumentError("Use ntheta≥8 and nr≥5"))
    isfinite(min_peak_contrast_fraction) && 0 <= min_peak_contrast_fraction < 1 || throw(ArgumentError("min_peak_contrast_fraction must lie in [0,1)"))
    1 <= m_max < ntheta/2 || throw(ArgumentError("m_max must be positive and below the angular Nyquist mode"))
    rlo, rhi = max(0.0, jet_radius_m - jet_halfwidth_m), Float64(jet_radius_m + jet_halfwidth_m)
    # Samples must fit inside one copy of the grid; the sampler still handles
    # interpolation across the last-to-first cell without a duplicated endpoint.
    extent = min(-grid.x0, grid.x0 + grid.nx*grid.dx,
                 -grid.y0, grid.y0 + grid.ny*grid.dy)
    rhi < extent || throw(ArgumentError("Search annulus must lie inside the Cartesian domain about (0,0)"))
    radii = collect(range(rlo, rhi; length=nr))
    theta = collect(0:ntheta-1) .* (2π/ntheta)
    ut, ur, pv = fill(NaN,nr,ntheta), fill(NaN,nr,ntheta), fill(NaN,nr,ntheta)
    for k in eachindex(theta)
        st, ct = sincos(theta[k])
        for j in eachindex(radii)
            xx, yy = radii[j]*ct, radii[j]*st
            uu, vv = _sample(u,xx,yy,grid), _sample(v,xx,yy,grid)
            ut[j,k], ur[j,k] = -st*uu + ct*vv, ct*uu + st*vv
            pv[j,k] = _sample(totalpv,xx,yy,grid)
        end
    end
    mean_utheta, mean_ur, mean_pv, eddy_uv = fill(NaN,nr), fill(NaN,nr), fill(NaN,nr), fill(NaN,nr)
    radial_coverage = zeros(nr)
    for j in eachindex(radii)
        good = findall(k -> isfinite(ut[j,k]) && isfinite(ur[j,k]), eachindex(theta))
        radial_coverage[j] = length(good) / ntheta
        if !isempty(good)
            tv, rv = ut[j,good], ur[j,good]
            mean_utheta[j], mean_ur[j] = mean(tv), mean(rv)
            eddy_uv[j] = mean((tv .- mean_utheta[j]) .* (rv .- mean_ur[j]))
        end
        goodpv = filter(isfinite, pv[j,:])
        mean_pv[j] = _finite_mean(goodpv)
    end
    jet_radius_theta, peak_speeds = fill(NaN,ntheta), fill(NaN,ntheta)
    dr = radii[2]-radii[1]
    for k in eachindex(theta)
        profile = view(ut,:,k)
        all(isfinite, profile) || continue
        scale = maximum(abs, profile)
        tol = max(1e-10, min_peak_contrast_fraction*scale)
        maximum(profile)-minimum(profile) > tol || continue
        j = argmax(profile)
        1 < j < nr || continue
        profile[j] > tol || continue
        fm, fc, fp = profile[j-1], profile[j], profile[j+1]
        curvature = fm - 2fc + fp
        curvature < -max(100eps(Float64)*scale, tol*1e-3) || continue
        offset = clamp(0.5*(fm-fp)/curvature, -1.0, 1.0)
        jet_radius_theta[k] = radii[j] + dr*offset
        peak_speeds[k] = fc - 0.25*(fm-fp)*offset
    end
    good = findall(isfinite, jet_radius_theta)
    coverage = length(good)/ntheta
    meanjet, mean_peak_speed = _finite_mean(jet_radius_theta[good]), _finite_mean(peak_speeds[good])
    jet_detected = length(good) == ntheta
    mode_coefficients = fill(complex(NaN,NaN), m_max)
    if jet_detected
        displacement = jet_radius_theta .- meanjet
        for m in 1:m_max
            mode_coefficients[m] = sum(displacement .* cis.(-m .* theta)) / ntheta
        end
    end
    mode_amplitudes = 2 .* abs.(mode_coefficients)
    return (; radii, theta, mean_utheta, mean_ur, mean_pv, eddy_uv,
            jet_radius_theta, mode_coefficients, mode_amplitudes,
            meanjet, mean_peak_speed, coverage, radial_coverage, jet_detected)
end

"""
    load_wind_profile(path; required_radius_range=nothing)

Read a two-column CSV with header `radius_m,u_theta_m_s`. Blank lines and
whole-line `#` comments are allowed. No quoted/embedded-comma fields or implicit
unit conversions are accepted. Require at least two finite samples, nonnegative
strictly increasing radii, and optional coverage of `(rmin,rmax)` in meters.
The returned named tuple contains the two vectors and their available interval.
This loader does not infer a reference frame, latitude projection or units.
"""
function load_wind_profile(path::AbstractString; required_radius_range=nothing)
    radius_m, u_theta_m_s = Float64[], Float64[]
    header_seen = false
    for (lineno, raw) in enumerate(eachline(path))
        line = strip(replace(raw, '\ufeff' => ""))
        (isempty(line) || startswith(line,"#")) && continue
        cells = strip.(split(line, ','))
        if !header_seen
            cells == ["radius_m", "u_theta_m_s"] || throw(ArgumentError("CSV header must be radius_m,u_theta_m_s (line $lineno)"))
            header_seen = true
            continue
        end
        length(cells) == 2 || throw(ArgumentError("Expected two CSV fields at line $lineno"))
        r, speed = tryparse(Float64,cells[1]), tryparse(Float64,cells[2])
        (r !== nothing && speed !== nothing && isfinite(r) && isfinite(speed)) || throw(ArgumentError("Nonfinite or invalid number at CSV line $lineno"))
        r >= 0 || throw(ArgumentError("Radius must be nonnegative at CSV line $lineno"))
        (isempty(radius_m) || r > last(radius_m)) || throw(ArgumentError("CSV radii must increase strictly (line $lineno)"))
        push!(radius_m,r)
        push!(u_theta_m_s,speed)
    end
    header_seen && length(radius_m) >= 2 || throw(ArgumentError("Wind CSV must contain a header and at least two samples"))
    if required_radius_range !== nothing
        length(required_radius_range) == 2 || throw(ArgumentError("required_radius_range must contain (rmin,rmax)"))
        lo, hi = required_radius_range
        all(isfinite, (lo,hi)) && 0 <= lo <= hi || throw(ArgumentError("Required radius interval must be finite, nonnegative and ordered"))
        first(radius_m) <= lo && last(radius_m) >= hi || throw(ArgumentError("Wind CSV does not cover requested radius interval [$lo,$hi] m; available [$(first(radius_m)),$(last(radius_m))] m"))
    end
    return (; radius_m, u_theta_m_s, radius_min_m=first(radius_m), radius_max_m=last(radius_m))
end

end # module
