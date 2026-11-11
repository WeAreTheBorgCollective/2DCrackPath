# CrackMesh2 / CrackStraight run

Mesh: 98 nodes, 160 triangles. Fixed comparison mask: 143175 pixels.

Status: no-improving-geometry. Accepted updates: 8.

| Stage | RMS [GL] | Mean absolute [GL] | Max absolute [GL] |
|---|---:|---:|---:|
| Regular DIC | 144.827132 | 93.451161 | 2293.143066 |
| XDIC step 0 | 129.618216 | 85.856720 | 2816.263989 |
| XDIC step 1 | 127.331410 | 85.282340 | 2823.915557 |
| XDIC step 2 | 123.643574 | 84.501306 | 1938.931265 |
| XDIC step 3 | 123.386980 | 84.429873 | 2121.574067 |
| XDIC step 4 | 123.298181 | 84.395277 | 2121.574067 |
| XDIC step 5 | 123.274970 | 84.392833 | 1938.931265 |
| XDIC step 6 | 123.253524 | 84.387239 | 2121.574067 |
| XDIC step 7 | 123.236830 | 84.379602 | 2121.574067 |
| XDIC step 8 | 123.195280 | 84.366765 | 2121.574067 |

ResidualComparison.png uses absolute original-image residuals with shared limits [0, 2823.915557] GL; high values are black. No clipping or per-panel normalization.

CrackPath_XX.csv contains reference coordinates and both deformed-image traces. Step_XX.mat saves the three level sets, coefficients and signed residual; Evolution.mat also saves rejected trials. CrackEvolution.png shows paths and RMS.

The mesh is translated by [imageWidth/2,imageHeight/2-90] as in example2Dscreening. The supplied crack is already in reference image coordinates and is not translated. Regular DIC uses sigma 2 then original images. XDIC evolution uses original images with a 5x5 residual screening window.


Validation run: 2026-11-11, MATLAB R2024b. The straight path length was 347 px; the final path length is 357.422642 px. The final solve uses 15 enriched nodes. Residual RMS falls by 4.96% from straight-crack XDIC and 14.94% from regular DIC on the same mask. No further improving geometry was found after eight accepted updates; this is not proof of a globally optimal or physically exact crack.

The sigma-2 paper PF run converged at iteration 11. All 10 PF tests passed, including the scale schedule and original-image residual checks.
