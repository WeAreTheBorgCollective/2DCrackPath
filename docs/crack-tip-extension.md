# Endpoint extension and final XDIC displacement

extractResidualCrackPath retains valley-Otsu thresholding, skeletonization and
graph selection. It then extends both ends of the selected ordered path.
The outward direction is the chord from a point 10 pixels inward along the
path to the tip (tipDirectionLength is editable). The two directions are
computed before either extension is appended.

For each tip, candidates are foreground pixel centers in its 8-connected
binary component, strictly ahead of the tip and at most 1 pixel from the
outward line. The largest forward projection selects the endpoint; ties
prefer the closest pixel to the line. The actual selected pixel is connected
to the path by a straight segment sampled at intervals no longer than 1 px.
The candidate rule can cross holes or concavities in the same component; it
does not connect to disconnected foreground islands. A tip outside the
binary region or without forward candidates is left unchanged.

CrackPathInit.mat now includes the pre-extension skeletonPathXY, binaryCrack,
crackSkeleton and tipExtension metadata. The usual crackPathXY/reference/
deformed fields and raster mask contain the extended result. graphInfo
continues to describe the original skeleton graph; tipExtension gives both
extension lengths and the final total path length. Regenerate extraction
before rerunning evolution to use these new initial tips.

runXDICCrackMesh2 is now a postprocessor of XDIC/Evolution/Evolution.mat.
Its historical filename is retained, but it no longer loads CrackMesh2,
starts from CrackStraight, correlates images or evolves a new crack.
Set xdicEvolutionFolder to read another evolution directory.
The field is evaluated as saved final pixel basis times saved final
coefficients, including Heaviside enrichments. Geometry and enrichment
ordering must match the saved basis. No current initialization or image
pair is consulted, so changing the active screening case does not mix runs.

FinalDisplacement.mat contains U (rows x columns x 2), ux, uy, magnitude,
crackPathReferenceXY and source information. Values outside the prepared
mesh pixels are NaN. FinalDisplacement.png shows the three displacement
maps with the last accepted crack overlaid in reference coordinates.
The black line is the reference crack geometry, not a deformed crack face.

Validation: horizontal and diagonal tips, reversed ordering, the 1-pixel
boundary, disconnected distractors, and already-complete paths passed.
Synthetic enriched-field reconstruction and mismatched-geometry rejection
passed. Reconstruction from the existing Evolution.mat matched its saved
kernel displacement within 1e-10 px. On the current original-image residual,
endpoints moved from (236,462)/(517,467) to (224,462)/(520,468): extensions
12.000/3.162 px, total length 315.630 px instead of 300.468 px. This validates
the geometric rule, not physical crack-tip accuracy. The trial did not replace
the active CrackPathInit.mat or rerun evolution.
