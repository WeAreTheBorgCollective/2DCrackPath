# Normal-line residual refinement

extractResidualCrackPath now refines the ordered principal Euclidean-graph
path after the existing binary-mask endpoint extension. Each ordered path
node (including added endpoint samples) gets an orthogonal search line from
its adjacent-node tangent. Endpoint tangents are one-sided. Directions are
computed from the unchanged input path; all node movements are simultaneous.
Discarded branches of the original skeleton graph are not crack-path nodes.

Edit normalSearch.radius (default 30 px) and normalSearch.density (default
1 sample/px). Offsets are (-floor(radius*density):floor(radius*density))/density.
Thus radius=30,density=1 gives exactly 61 samples, and density=2 gives 121.
Fractional products stop at the last sample within the requested radius.
Radius=0 disables movement.

The score is absolute RAW Residual/R.mat, not the screened rho used for
segmentation. Values are bilinearly interpolated over finite ROI neighbors.
Outside-ROI and outside-mesh query points are excluded. Each node selects
the actual largest valid sampled value. Ties prefer the nearest sample,
keeping the original node on flat profiles. Equal-distance ties use offset
order. No smoothing or smaller-peak substitution is applied.

The selected points replace crackPathXY and its reference/deformed mappings.
refinedCrackGraph is the new ordered chain with recomputed Euclidean weights
and original coordinates/normal offsets on its nodes. crackGraph and graphInfo
retain the original skeleton graph diagnostics. preNormalPathXY and
normalSearchInfo save the original path, sample offsets/scores, selected
offsets, before/after residuals and refined length. The path mask and plots
use the refined coordinates.

Independent maxima can create crossings, collapsed edges or links crossing
the ROI. The exact requested maxima are retained; validForLevelSets records
the domain/self-intersection check and extraction warns if it fails. Inspect
these trials before XDIC. This is an initialization experiment, not a guarantee
of a valid crack geometry. Clear cached rho when changing images/cases.

Validation: synthetic tests passed for horizontal/vertical normals, endpoint
directions, exactly 61 samples at radius 30/density 1, increased density,
radius limits, absolute residual polarity, flat-profile ties and ROI exclusion.
The complete extraction/mapping/graph-save/plot workflow ran in an isolated
folder. On the saved original-image path, 226 of 282 nodes moved; the raw
independent-maxima path failed the geometry-validity check. The full extraction
trial also flagged invalid geometry. These outputs were inspected and the
active initialization was preserved during validation.
