# ROI screening and crack-face length updates

screenResidual2D now evaluates every finite ROI center, including image and
ROI boundaries. It sums absolute finite residuals in the clipped rectangular
window and divides by the count of those samples. Outside pixels contribute
neither a value nor a count, and outside-ROI centers stay NaN. Pass an optional
sixth logical roiMask to exclude finite values outside an explicit ROI;
existing callers use NaN/Inf exclusion. Clear cached rho before extraction.
This removes the half-window boundary loss; it cannot recover pixels absent
from the DIC kernel's measurement domain.

evolveXDIC now alternates its shape proposal/solve with a face-growth proposal.
Growth is attempted even if the shape proposal stalls. The two physical
deformed-image traces are calculated from the current enriched solution.
At each endpoint, outward directions use chords spanning 10 reference pixels
along the path. A unique intersection must be strictly forward along BOTH
rays. Parallel rays, coincident tips and intersections behind either ray
produce no growth.

The intersection is in the deformed frame. Its increment from the midpoint
of the two face tips is mapped back with the inverse local deformation
gradient of the mean trace. This gradient comes from the endpoint triangle's
affine nodal standard and enriched displacement, and is extrapolated locally.
The midpoint here defines a geometric mapping convention, not a unique
physical displacement on the crack. Singular/inverted local mappings are
rejected. Reference increments must point outward and are capped at 10 px.

All three level sets are rebuilt for a growth proposal. Enrichment is rebuilt,
coefficients are transferred by node ID and XDIC is solved again. The existing
fixed-mask residual-improvement rule still applies. Out-of-domain or crossing
paths are rejected; an intersection alone does not force acceptance.
Shape optimization resumes from accepted growth on the next outer iteration.
The outer iteration limit bounds both stages; each accepted stage has its own
consecutively numbered MAT/PNG snapshot and history entry.

Configure options.faceGrowth in runXDICCrackEvolution:
enabled=true, maxExtension=10, tangentLength=10.
Set enabled=false for shape-only behavior. Optional parallelTolerance=1e-6
is the minimum absolute determinant of the normalized ray directions.
result.attempts records tip intersection coordinates, ray parameters,
uncapped/capped reference growth and solve/acceptance reasons.

This is the requested geometric heuristic, not a fracture propagation law or
a claim of physically exact tip position. Current Heaviside traces may
already meet at the geometric tip; that gives zero forward growth.

Validation (2026-11-17): 13 screening/Rethore regression tests and the XDIC
level-set checks passed. Ray tests cover forward/parallel/coincident cases,
the 10 px cap, reversed ordering and affine frame conversion. The original
dataset now has rho at all 407035 finite residual pixels. A two-iteration
trial accepted shape, growth, shape and growth in that order: RMS values
114.8572 -> 114.8535 -> 114.6952 -> 114.6297 -> 114.6223 GL. Both tips grew
10 px at each accepted growth stage. Saved paths, basis reconstruction and
cap checks passed. This bounded trial is not a claim of convergence or
physical tip accuracy; it leaves the active evolution output unchanged.
