# Two-dimensional triangular XFEM DIC

Run `buildXDICKernels` once with a configured MATLAB C++ compiler
(`mex -setup C++`), then `example2DXDIC`. Correli must be on the MATLAB path.
The example reads the images, mesh and regular initialization used by
example2Dscreening and CrackPaths/CrackPathInit.mat. It saves the comparison
to XDIC/Comparison.mat and exports residual, displacement and enrichment plots.
Run `addpath tests; testXDIC` for validation.

## Formulation

Adaptation of the shifted-Heaviside formulation in Réthoré,
Extended digital image correlation with crack shape optimization,
Section 2.3, equations 15–18, to linear triangles and pixel quadrature:

u(X) = sum_i N_i(X) a_i + sum_j N_j(X) [H(X)-H(X_j)] b_j.

Only nodes whose support intersects the finite crack are candidates.
Every node belonging to an element containing either endpoint is excluded.
No square-root crack-tip functions are included. Both sides must have
nonzero sampled support; unresolved candidates are reported in preparation.
This targets a single simple open crack on a conforming T3 mesh.
Branching, self-intersecting cracks and overlapping/degenerate meshes are
unsupported. Endpoints outside the mesh allow a through-crack.
H is 0/1 from the nearest oriented crack segment; reversing the path
reverses enriched coefficient signs without changing the approximation space
away from exact on-crack ties (assigned H=1).

All coordinates are reference-frame, one-based image pixels (x=column,
y=row). Do not pass the mapped deformed-frame crack.
Residual is f(X)-g(X+u(X)); no damage is applied.

## Interfaces

- prepareXDIC(Mesh,imageSize,crackPathXY) wraps xdicPrepMex.
  Mesh.pos has XY or planar XYZ coordinates; Mesh.con has [0 n1 n2 n3].
  Direct MEX inputs are full doubles: Nx2 nodes, Ex3 connectivity,
  [height width], and crack Nx2 (or empty).
- prep.basis is sparse pixels-by-(nodes+enrichedNodes).
  prep.pixelIndex maps rows to MATLAB image indices.
  Standard and enriched basis values, node lists, tip supports and cut
  elements are available for inspection. Each pixel is integrated once,
  including shared triangle edges. Outer image-border pixels are excluded.
- xdicCorrMex(prep,f,g,a,interpolation,assemble) accepts single/double images
  and double coefficients interleaved ux,uy for standard nodes followed by
  enriched nodes. Interpolation is linear or cubic (four-point Lagrange).
  Returns H, b, R, RMS r, standard nodal residual sums, U, validity and counts.
  With assemble=false it evaluates residual/displacement only.
- correlateXDIC controls the MATLAB solve, diagonal scaling, damping,
  backtracking and pixel-increment convergence. It returns coefficients,
  standard/enriched amplitudes and the full image displacement.
  Empty crack gives ordinary T3 FE with zero enriched DOFs.

The tangent uses fixed reference central-difference gradients, matching
Correli's modified Gauss–Newton kernel: J=[B*f_x,B*f_y], H=J'J, b=J'R
(with interleaved columns). All mixed standard/enriched blocks are assembled.
Invalid interpolation pixels are skipped. Line search preserves initial valid
coverage and requires nonincreasing residual there. Unobservable DOFs are
held fixed. No mechanics, brightness correction, 3D or tip enrichment is used.
Preparation and assembly are standalone C++ MEX; linear algebra/control remain
MATLAB. Existing Correli files and binaries are not changed.

## Validation on MATLAB R2024b / MinGW64

Tests cover partition of unity, affine reproduction, tip exclusion, reversed
crack orientation, explicit Jacobian comparison for every Hessian/RHS block,
and synthetic 0.6 px jump recovery. Synthetic RMS: ordinary FE 0.312355,
XFEM 0.0021048; displacement RMS error 0.000504 px.
Empty-crack assembly against installed Correli on a noninteger triangle:
relative Hessian error 5.38e-16, RHS error 2.93e-14,
maximum residual discrepancy 3.13e-13.

Project comparison with identical initialization and common valid pixels:
Correli RMS 122.639, new FE RMS 122.639, XFEM RMS 117.175 (4.46% reduction).
11 enriched nodes; both new solves reported convergence.
This validates implementation behavior, not accuracy of the supplied crack.
The fixed crack and omitted tip enrichment limit fitting near endpoints.
Compiled MEX files and generated MAT/PNG results are local build artifacts.
