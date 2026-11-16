# Extracted-crack initialization, original images

Run extractResidualCrackPath first, then runPaperPhaseFieldDIC.

The PF script loads CrackPaths/CrackPathInit.mat and uses its explicitly
mapped crackPathDeformedXY. The paper solver operates on the deformed image
grid, so the reference-frame crack is not used directly.

At each ROI pixel, the initial damage is
d0=exp(-sqrt(mu)*distanceToFiniteCrack/lc).
This is a diffuse distance-profile initial guess matching the homogeneous
straight-crack damage operator; it is not an exact curved-crack PDE solution.
Outside the measurement ROI the saved initial damage is zero. The existing
lc and mu control the seed width (default decay length 15.81 pixels).
The profile seeds the first weighted displacement solve and is then free
to change under the independent damage PDE. No crack constraint or damage
irreversibility has been added. Initial displacement remains zero.

No Gaussian image filter is used in this entry point, even if a stale
paperParameterOverrides specifies smoothing. Displacement and damage use
the residual g(x)-f(x-u(x)) from the original images. Existing affine
graylevel normalization, intensity inversion and the paper's block-based
Phi calculation remain; the raw squared residual term itself is not screened.

Every timestamped run saves InitialDamage.mat, the mapped source path,
profile metadata in Model.mat, Iteration_00.mat/png and subsequent iterations.
Missing deformed coordinates, mismatched image sizes and out-of-domain paths
produce errors instead of silently using the wrong frame or a zero seed.
This supersedes the earlier zero-damage and Gaussian initialization defaults.
