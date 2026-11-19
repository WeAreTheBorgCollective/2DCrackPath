# Fixed-field gradient refinement during extraction

The independent normal-line maximum search from PR #21 is removed.
Extraction again starts with the previous valley-Otsu, skeleton graph and
binary-mask endpoint-extension stages, then refines that curve by gradient
level-set evolution on a frozen regular-DIC score field.

In extractResidualCrackPath, choose fixedFieldSource='rho' (default) to use
the existing screened residual, or 'residual' for abs(Residual/R.mat).
The chosen field and ROI are copied once and never refreshed. A supplied rho
is not screened again. Clear cached rho when changing the DIC input pair.

The shape updater uses normal velocity proportional to grad(field).normal,
the existing Godunov equation, tangential speed smoothing and backtracking.
Defaults: 30 iterations, gain 0.1, maximum normal speed 0.5 px per virtual time,
smoothing span 3 px, motion tolerance 0.02 px. Only valid, non-crossing,
mesh-contained curves that increase the mean frozen score are accepted.
There is no correlation re-solve, global-RMS veto, face-intersection growth
or element-tip snapping in this stage. Endpoints follow normal shape motion.
This is a local gradient fit, not a guaranteed global best path.

CrackPathInit.mat contains the refined reference/deformed path and graph.
The original graph diagnostics are retained; preRefinementPathXY records the
starting curve. FixedFieldRefinement.mat stores the frozen field, final level
sets, every accepted path, score history, options and attempt diagnostics.
FixedFieldRefinement.png shows path progression and the accepted score.
No changes are made to runXDICCrackEvolution or its existing solver behavior.

The propagated zero contour is resampled as in the existing XDIC shape
updater, so the refined graph can have a different vertex count. Synthetic
ridge and flat-field tests passed. A five-step saved-data trial increased
mean fixed score from 254.779 to 256.450, with exact frozen-field equality
and strictly increasing accepted scores verified. These checks establish
local ascent and data consistency, not global optimality.

An isolated three-iteration end-to-end extraction also passed: score
251.798 -> 252.696, exact frozen-field equality and valid final geometry.
Graph persistence, coordinate mapping and progression plots completed.
The production default remains 30 iterations; active initialization files
were not overwritten during validation.
