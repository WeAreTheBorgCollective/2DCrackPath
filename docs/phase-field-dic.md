# Phase-field DIC and crack-level-set evolution

## Run order

Configure Correli, then run:

```matlab
example2Dscreening
screenResidualImage
extractResidualCrackPath
runPhaseFieldDIC
```

The initial screening window is 23 x 23. Clear an existing workspace `rho`
before standalone extraction if you want it recomputed with the new default.

`runPhaseFieldDIC` is the sole PF-DIC driver. The redundant
`refineCrackPathRethore.m` script has been removed. The driver starts from
`CrackPaths/CrackPathInit.mat` and `InitialValues/DisplacementInit.mat`.

## Order of operations

1. Construct the crack-surface level set phi and two tip level sets from the
   extracted polyline, then construct its damage field.
2. Call PF-DIC with that damage and the current displacement initialization.
3. With the accepted displacement, evaluate
   `R(X) = f(X) - g(X + uPF(X))` using the original reference and deformed
   images and the native Correli interpolation kernel. **No damage is applied
   to either image during this residual evaluation.** The dedicated
   `evaluateUndamagedResidual` helper has no damage argument.
4. Screen this residual, evolve the level sets using the normal component of
   its screened gradient, then reconstruct all three distance fields and damage.
5. Save the state and repeat. A final PF-DIC solve after the last shape update
   synchronizes the saved displacement with the final damage.

Displacement proposals retain the energy acceptance check described below.
A rejected proposal keeps the previous accepted displacement; snapshots and
diagnostics do not claim that rejected candidates are PF solutions.

## Level sets and damage

`crackLevelSets2D` constructs phi as signed distance to the oriented polyline,
extended beyond both tips by tangent rays. Positive phi is left of the directed
curve in image coordinates. The tip fields are
`psiStart = -(X-Pstart).tStart` and `psiEnd = (X-Pend).tEnd`.
Their zero contours intersect the crack at its tips and are locally orthogonal
to phi. Exact unsigned distance to the finite polyline controls damage, so the
tangent extensions do not create damage beyond the tips.

These are the two-tip adaptation of Rethore, Hild and Roux,
*Extended digital image correlation with crack shape optimization*,
DOI 10.1002/nme.2070, Eqs. 20-22. Local tip planes can intersect other portions
of folded cracks; finite-polyline distance remains authoritative.

Damage inside the fixed measurement ROI is the geometry-constrained profile

```text
d = exp[-sqrt(mu)/lc * (sqrt(distance^2 + epsilon^2) - epsilon)].
```

Default lc=5, mu=0.1 and epsilon=0.1 px. Outside-ROI values are preserved
(default zero). This is a distance-based profile, not a solution of the
independent residual-driven damage PDE in the PhaseFieldDIC paper.

## Screened-gradient motion

Configure the independent level-set window in `runPhaseFieldDIC`:

```matlab
levelSetOptions = struct( ...
    'screeningWindowSize',[23 23], ...
    'velocityMethod','gradient', ...
    'gradientGain',1, ...
    'maxNormalStep',1, ...
    'duration',1, ...
    'smoothLength',3);
```

Odd rectangular windows such as [7 15] are supported. [1 1] uses the absolute
raw residual. Each step recomputes rho = mean(abs(R)) over the chosen centered
window. Only complete finite windows inside the ROI are eligible.

Normal velocity is `F = clip(gradientGain * grad(rho).normal, +/-maxNormalStep)`,
smoothed along arclength. Gain has units px^2/GL per virtual time, the speed cap
has units px per virtual time, and duration controls virtual time. Unsupported
gradient stencils are not used, and smoothing cannot activate unsupported nodes.
The Godunov solver integrates `phi_t + F*|grad(phi)| = 0` with CFL substeps.
Backtracking accepts only an increase of the whole-path screening score and
valid connected geometry. Unsupported window centers contribute zero to that
score, with no raw-residual fallback.

This local gradient-ascent velocity is an explicit adaptation, **not**
Rethore's Eqs. 35-36 normal-line maximum-search rule. The latter remains available
as `velocityMethod='maximum'` for comparison. Both use the Godunov integrator
of Eqs. 26-32. Exact distance reconstruction from the finite zero contour
replaces the paper's PDE reinitialization. Both tip fields are rebuilt after
each accepted update; no independent tangential crack-growth law is imposed.

PF energy is diagnostic for geometry; screening-score ascent need not decrease it.

## Correlator limitation and acceptance

`correlatePhaseField` uses the existing `correlate` on images multiplied by
(1-d), as requested in the original workflow. This is an image-masking surrogate:

```text
masked residual = (1-d(X))*f(X) - (1-d(X+u))*g(X+u)
desired weighted residual = (1-d(X)) * [f(X)-g(X+u)]
```

These differ for nonzero displacement and nonuniform damage. Correct workflow
ordering and raw-residual evaluation do not remove this limitation. An exact
weighted PF-DIC solve would require weighting inside the correlator and its
linearization. BCC and mechanical regularization are disabled.

Proposed displacements are reevaluated on the original images and accepted only
if they preserve ROI coverage and do not increase the mechanical-free energy:

```text
E = 1/2 sum_ROI [(1-d)^2*(R/grayScale)^2 + (mu/lc)*d^2]
  + lc/2 sum_internal_bonds (d_j-d_i)^2.
```

The fixed gray-level scale comes from the reference image. Only bonds with
both pixels in the ROI contribute. Damage is geometry-constrained; the code
does not implement a mechanical driving term, independent damage PDE, crack
nucleation, or irreversibility.

## Outputs and visualization

Each run creates `Damage/Evolution_<timestamp>` and preserves the prior named
damage/result MAT files there. It saves `Step_00` (initial), one snapshot per
shape step, and a final PF-solve snapshot. Each MAT contains all three level
sets, damage, accepted displacement, unmasked residual, screened residual,
options and a stage label. Each PNG shows the raw residual, screening, phi
with tip zero contours, and damage. After a shape update the residual belongs
to the preceding PF solve; the final snapshot follows the final PF solve.

Named final outputs are `Damage/Damage.mat`, `Damage/PhaseFieldResult.mat`,
and `InitialValues/DisplacementInitPF.mat`. The driver automatically runs
`comparePhaseFieldDIC` to refresh the initial/final displacement and original-
image residual comparisons. Histories record geometry and DIC acceptance
separately, including the final DIC solve.

## Validation

```matlab
results = runtests('tests');
assertSuccess(results);
```

Tests cover distance geometry, damage/energy derivatives, normal-gradient
direction and speed cap, independent screening windows, ROI holes, stationary
profiles, the Godunov propagation sign, original-image kernel arguments,
displacement packing and rejected residual coverage loss.

On the repository data, the consolidated driver completed three geometry steps
and a final PF solve. The first gradient update increased the screening score
from 259.59807 to 262.01905 GL; the next two were rejected. All four DIC
proposals were rejected, so the final displacement and original-image residual
were unchanged (RMSE 122.64 GL). Five MAT/PNG snapshots were generated.
All 29 MATLAB tests passed. These results expose the remaining image-masking
surrogate limitation rather than demonstrating a converged weighted PF solve.

## Crack coordinates and deformed-image overlays

The residual is evaluated at reference pixels X:
`R(X)=f(X)-g(X+u(X))`. Extracting a crack from this image therefore already
gives reference coordinates. Moving these coordinates a second time for the
damage/level-set computation would mix frames.

For display on the original deformed image, map every crack node forward:
`x=X+u(X)`. The displacement is interpolated using barycentric shape functions
of the actual Correli T3 connectivity, preserving mesh holes; unsupported points
raise an error rather than being extrapolated. Subtracting u at X is not the
appropriate mapping for this residual convention.

`CrackPathInit.mat` now stores `crackPathReferenceXY`,
`crackPathDeformedXY`, and `crackPathDisplacementXY`, together with frame and
mapping metadata. The compatible `crackPathXY` field and mask remain in the
reference frame for PF-DIC. Extraction uses
`InitialValues/DisplacementInit.mat`, the displacement associated with the
initial residual, and saves `CrackPaths/CrackOnDeformed.png` and
`CrackPaths/CrackExtraction.png`.

Each evolution snapshot also stores both crack coordinate sets and overlays
the mapped path on the original grayscale deformed image. The mapping uses
that snapshot's accepted displacement a, rather than always using the initial
or final field. The raw-residual, screening, level-set and damage panels remain
in reference coordinates. No damage is applied to the displayed deformed image.
