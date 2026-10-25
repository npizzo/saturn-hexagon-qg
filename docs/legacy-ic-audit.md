# Audit of the supplied legacy initializer

The supplied MATLAB text contains a useful Ptolemaic PV formula. Its **unshifted
`q0` matches the independently differentiated planar map**, but the meaning of
the saved `q.txt` in the original spherical solver remains unverified. This audit
does **not block the existing annular-jet campaign** and does not certify a
complete exact, stable, or correctly drifting Ptolemaic QG initial state.

The attachment SHA-256 is
`cd2afd3bb7c391b512ba341fc8f8ee751412f5f07d49c1ebd019501f655c2b47`.
The local audit record includes the source path,
script hashes, all 16 comparison points, active parameters, and unresolved
conventions. MATLAB and the legacy solver were not executed; no production
initializers or model files were changed.

## What the active code does

The attachment sets `n=2`, `eps=0.05`, and `a0=0.3` (lines 16–20). The map has
polygonal symmetry `m=abs(n+1)=3`; a hexagon uses `n=5`. It samples
`b=-2:0.0025:0` and constructs a grid expected to be 201×201. Its `omega=5` is
planetary rotation because `gamma=2*omega`; the assigned `c=0` is unused.

After building the map and `q0`, the code resets `eps=0`, but the calculation of
the axisymmetric comparison `q1` is commented out. Thus `q1` stays zero and no
axisymmetric subtraction occurs (lines 45–56).

The active pipeline interpolates `q0` with nearest neighbors, subtracts the
gridded minimum, replaces NaN/Inf with zero, and saves `Q=q2` (lines 68–95, 146).
The spherical Coriolis field `f` and mask `q3` are computed but not used in `Q`.
This is a shifted, truncated, interpolated PV-shaped patch. A minimum subtraction
is not automatically harmless when the result is passed to a vorticity inversion.

MATLAB's `griddata(...,'nearest')` returns NaN outside the sample convex hull;
the script subsequently sets those values to zero. SciPy nearest-neighbor
extrapolation alone would not reproduce that exterior. MATLAB `min` omits NaN
by default. [griddata documentation](https://www.mathworks.com/help/matlab/ref/griddata.html),
[min documentation](https://www.mathworks.com/help/matlab/ref/double.min.html).

## Formula verification

Define `p=exp(b)`, `h=eps*exp(n*b)` and `D=p^2-n^2*h^2`. The attachment's expression
is

```text
q0 = gamma * [1 - A0^2 * ((n+2)*p^4 - n^3*(1+2*n)*h^4)/(n*D)].
```

It agrees with total planar barotropic PV `qBT=laplacian(psi)+f`, under
`u=-psi_y`, `v=psi_x`, `f=gamma-gamma*r^2`, and

```text
W(b) = -gamma*A0^2*(p^2-n*h^2)/(2*n),
theta = a-W(b)*t.
```

The additive wave-frequency constant is zero in this matching construction;
it is distinct from the script's planetary `omega`. The maximum absolute
difference at 16 `(a,b)` points is **1.78×10⁻¹⁵**. The comparison uses the
map-derived analytic curl in `audit_ptolemaic_map.py`, whose Cartesian
finite-difference verification is recorded separately in
`validation/ptolemaic-map-audit.json`. The supplied formula therefore supports
the audited opposite-to-printed frequency sign under these explicit conventions.

## What remains unresolved

The named folder
`Rick.Sphere.Code/vortex.waves`
returned macOS **Operation not permitted** during targeted inspection. Neither
its reader nor its inversion was inspected. Filenames and hemisphere comments
cannot establish whether `q.txt` is relative vorticity, absolute/total vorticity,
or a solver-specific anomaly. Inertial relative vorticity can represent absolute
vorticity in a rotating description, so the reference frame also matters.

The plotting map uses northern `z=(1-rho^2)/(1+rho^2)`. The computed but unused
`f=2*Omega*(1-rho^2)/(1+rho^2)` expands as `2*Omega-4*Omega*rho^2`, while the
verified planar formula uses `2*Omega-2*Omega*r^2`. Their coordinates therefore
need an explicit scale/chart conversion before claiming equivalence. Under the
plotted stereographic map, `a0=0.3` corresponds to latitude 56.6°, with the
outer contour spanning approximately 55.0–58.2°; these are not calibrated Saturn
hexagon parameters.

Before reusing this as a physical Ptolemaic initial state, establish the reader's
PV and mean conventions, hemisphere signs, coordinate/time scaling and array
ordering, and the complete outer patch. Then check the resulting velocity and
equilibrium residual. The already-defined annular-jet experiments can proceed
independently of that historical reconstruction.

## Reproduce

From the project directory, using Python 3.8 or later:

```sh
python scripts/audit_legacy_ic.py \
  '/path/to/legacy_initializer.m' \
  --output validation/legacy-ic-audit.json
```

The script reads the attachment as text, rejects an unrecognized `q0` expression,
and does not evaluate MATLAB or retry access to the legacy solver directory.
