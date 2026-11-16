# Sauvola residual segmentation trial

Run extractResidualCrackPathSauvola after example2Dscreening and
extractResidualCrackPath. It reads Residual/R.mat; it does not filter the images.
It compares valley-Otsu on the 23x23 screened residual, Sauvola on raw abs(R),
and Sauvola on the same 23x23 screened residual. The existing
CrackPaths/CrackPathInit.mat is not overwritten. Each trial writes a compatible
CrackPathInit.mat inside a timestamped CrackPaths/Sauvola directory when a
usable path exists. LatestRun.mat identifies the directory.

Sauvola was designed for dark document features:
T=m*(1+k*(s/R-1)), with local mean m and population standard deviation s.
To handle bright crack residuals, I=1-min(score/scale,1) and foreground is I<T.
The equivalent threshold in graylevels is scale*(1-T), compared with the
clipped score. This is a polarity adaptation of Sauvola, not a DIC-specific
crack-length estimator.

Defaults: 51x51 threshold window, k=0.2, R=0.5 on normalized [0,1] intensity.
The intensity scale is the finite-score 99.5 percentile, stored as scaleGL;
the clipped count is saved. The 23x23 screening window is a separate setting.
Finite ROI samples alone supply statistics. At least 50% of the requested
window must be finite; centers outside the ROI never become foreground.
This partial-window rule avoids treating NaNs as zero residual or eroding
every boundary by half the threshold window. Empty/constant inputs return
empty masks. This still leaves boundary bias and parameter sensitivity.

Each method uses the existing thinning and graph extraction, with identical
graph settings: MaxGap=25 px, MinComponentDiameter=25 px, MaxVertexDensity=2.
The shorter MaxGap differs from the original extraction script's 100 px.
Metrics.csv includes endpoints, principal path length, remaining components,
excluded nodes and the inferred bridge length actually on that path.
FullComparison.png shows all data; Comparison.png is a display-only zoom
around the previous crack. That previous crack is not used for segmentation
or path selection. Paths are mapped to the deformed image using the existing
regular DIC field; no XFEM solve is performed in this trial.

The raw residual MAT, masks, local threshold maps and parameter metadata are
preserved in Segmentation.mat. No physical length accuracy is inferred from
a longer, more connected or visually smoother result. A ground-truth tip
location or independent displacement-based check is still needed.

Reference: Sauvola & Pietikainen, Adaptive document image binarization,
Pattern Recognition 33 (2000), 225–236.
https://doi.org/10.1016/S0031-3203(99)00055-2

Validation: tests/testSauvolaResidual checks a bright synthetic ridge,
background rejection, graylevel-scale invariance, the local formula against
an independent calculation, missing pixels and constant inputs.

## Project-data trial (2026-11-16)

| Method | Foreground pixels | Retained path length [px] | Inferred length [px] | Reference endpoints |
|---|---:|---:|---:|---|
| Valley-Otsu, screening 23x23 | 8025 | 300.468 | 0 | (236,462) to (517,467) |
| Sauvola, raw absolute residual | 47451 | No retained path | 0 | None |
| Sauvola, screening 23x23 | 7584 | 359.379 | 21.840 | (236,462) to (563,452) |

The raw trial produced 15215 skeleton components, all removed by the graph
filters. This is not a measured zero crack length. The screened trial joins a
separate right-hand region and bends away from the main band. The extra
length therefore cannot be claimed as improved tip accuracy. These parameters
are exploratory, not calibrated. Existing production initialization remains
unchanged. The independently computed local-threshold and synthetic-ridge
tests passed.
