# CrackMesh2 and a straight initial crack

Run `runXDICCrackMesh2` with Correli on the MATLAB path and the XDIC MEX
kernels built. The script uses Mesh/CrackMesh2.m (98 nodes) and
CrackPaths/CrackStraight.mat, specifically crackPathReferenceXY.
The supplied line runs from (213,465) to (560,465) in reference image pixels.

The mesh is positioned with the existing example's translation
[imageWidth/2,imageHeight/2-90,0]. The crack is already in image coordinates
and is not translated. Its containment in the actual triangular mesh is
checked. Regular DIC is initialized from zero on the new mesh using Gaussian
sigma 2 followed by original images. No coefficients from the older mesh
are reused.

XDIC then solves on the straight crack and alternates residual-gradient
level-set proposals with XDIC re-solves. Settings are 12 maximum geometry
steps, 5x5 residual screening, gradient gain 0.1, maximum normal speed 1,
and tangential smoothing length 5 pixels. The three level sets and
rejection/backtracking logic are described in XDIC.md.

## Outputs

Each execution creates a timestamped directory below
XDIC/CrackMesh2Straight; LatestRun.mat points to it.

- Inputs.mat: positioned mesh, translation, supplied path and regular DIC.
- Step_XX.mat/png: accepted geometry, all three level sets, displacement,
  signed original-image residual and both mapped crack-face traces.
- CrackPath_XX.csv: reference x/y and deformed minus/plus x/y coordinates.
- ResidualMetrics.csv and Results.md: RMS, mean absolute and maximum
  absolute residuals in original graylevels, on one fixed common pixel mask.
- ResidualComparison.png: regular DIC, straight-crack XDIC and final XDIC.
  Every panel displays abs(R), using identical limits from zero to the
  maximum over all documented stages. White is zero; black is high.
  No percentile clipping, logarithm or image-intensity normalization.
- CrackEvolution.png: reference paths, deformed-image midpoint evolution,
  final one-sided crack-face traces and RMS history.
- Evolution.mat: final solution, paths, options and rejected-trial diagnostics.

The single-step zoom plots use their own full nonnegative residual range;
use ResidualComparison.png for direct comparisons with a shared scale.
Image dimensions and measurement domains differ from the CrackMesh1 example,
so residual metrics should not be compared across those ROIs as if identical.

## Gaussian initialization for paper PF-DIC

runPaperPhaseFieldDIC now defaults to displacementSmoothing=[2 0], matching
the Gaussian-then-original sequence in example2Dscreening. Both images are
filtered with imgaussfilt at sigma 2 pixels for the first displacement solve,
followed by an unfiltered refinement on the same grid. Subsequent outer
iterations use the originals. This is image smoothing continuation, not
spatial downsampling.

The damage PDE, its driving residual and final diagnostics use the
unfiltered images after refinement; displayed backgrounds are original
images. The sequence is saved in p and history(1).displacement.smoothingScales.
Override it with paperParameterOverrides.displacementSmoothing if needed.
The test suite verifies the [2 0] sequence, raw-residual evaluation and
unfiltered subsequent iterations.

Validated in MATLAB R2024b: all 10 paper PF tests passed. The full sigma-2
PF run converged after 11 outer iterations (weighted increment ratio 0.04608,
tolerance 0.05). This stopping criterion is not a physical damage calibration.

## Measured new-mesh result

See XDICMesh2Results.md for every accepted stage. On 143175 common pixels,
regular DIC RMS is 144.827132 GL, straight-crack XDIC is 129.618216 GL,
and the final result after eight accepted updates is 123.195280 GL.
The final reference path has length 357.422642 px (initially 347 px).
The run stops with no-improving-geometry. ResidualCrackZoom.png provides
a closer view with the same absolute graylevel limits as the ROI comparison.
