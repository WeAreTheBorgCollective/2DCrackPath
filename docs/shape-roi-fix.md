# Fix stalled residual-gradient shape updates at ROI borders

Reproduced with composite11 versus composite00 on mesh200 and mesh050:
the normal velocity was nonzero (reaching 1 px per virtual time unit), but
the before/after contour scores were -Inf and no shape was accepted.

The screening routine had been updated to use partial ROI windows, but the
shape scorer still required every bilinear neighbor of a path sample to be
inside the ROI. A valid fractional endpoint near the notch failed that test.
NaN finite-difference stencils also removed otherwise usable border gradients.
This was a boundary-sampling inconsistency, not a disabled level-set equation.

The update now:
- Computes central gradients inside the ROI and one-sided gradients at borders.
- Renormalizes bilinear interpolation over finite ROI neighbors, using nearest
  pixel membership for the query center, consistent with crackPathInDomain.
- Never introduces zeros or residual values from outside the ROI into sampling.
- Constrains subpixel contour excursions (at most 1 px) onto the actual mesh
  boundary before re-solving; a tiny inward numerical offset avoids edge
  roundoff. Larger excursions and invalid segments are rejected.
- Rechecks the screened score after mesh constraints and tip snapping.
- Distinguishes an accepted geometric proposal from global-RMS rejection in
  saved attempt reasons. Supported-point and contour-rejection counts are saved.

The gradient gain, Godunov propagation, screening window, face-growth rule and
global RMS acceptance criterion are unchanged. No MEX rebuild is needed for
this MATLAB-only change. Rerun runXDICCrackEvolution, then runXDICCrackMesh2.

Validation:
- ROI interpolation, one-sided gradients, fractional boundary endpoints,
  constant-field stationarity and mesh projection/rejection tests passed.
- Existing XDIC level-set tests and all five screening tests passed.
- Separate one-iteration shape-only trials (face growth disabled to isolate
  shape movement, maxIterations=40 for each XDIC solve) accepted and saved
  motion on BOTH composite cases. RMS fell from 3.9453207 to 3.9398527 GL on
  composite200 and from 3.4381947 to 3.3618547 GL on composite050.
- Proposed/saved motion was approximately 1.10 and 1.16 pixels respectively.
  These are one-step regression runs, not claims of full shape convergence.
  Active result folders were left unchanged; trial outputs are local.
