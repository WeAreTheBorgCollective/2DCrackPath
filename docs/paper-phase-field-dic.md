# Independent paper-based phase-field DIC

Run `runPaperPhaseFieldDIC` after configuring the project paths. It reads
the original images and only the bounding rectangle of the previous measurement
ROI. No extracted crack, damage profile, residual screening map, or prior
displacement initializes the solver. The old level-set experiment remains
separate in `runPhaseFieldDIC`.

Source: Zhu and Lu, *Phase-field digital image correlation for integrated
displacement and damage measurements*, arXiv:2606.16850v2, Sections 2.3-2.7.
The implementation below follows the staggered method, with numerical choices
and the registration safeguard stated explicitly.

## Equations and iteration order

The coordinate frame is the paper's deformed-image grid x, with
`r(u)=g(x)-f(x-u(x))`. Displacement u(x) maps back to reference coordinates.
It is NOT interchangeable with the older reference-grid Correli a(X).
Damage, its contours and skeleton are already in deformed-image coordinates
and are overlaid directly on the original grayscale g.

Initialize u=0 and d=0. For each outer iteration:

1. Hold d fixed and solve weighted image registration. At each Gauss-Newton
   step form `J=[diag(f_x(x-u))*N, diag(f_y(x-u))*N]`,
   `K=J' W J` and `Q=J' W r`, with `W=diag(mass*(1-d)^2)`.
   Solve for the displacement increment and update u. Images themselves are
   never multiplied by damage. This uses the shifted gradient in Eq.21 rather
   than the optional fixed-gradient approximation in Eq.22, and adopts the
   paper's approximation of dropping the Phi derivative in displacement.
   The printed Qu in Eq.31 omits the image-gradient factor; the implementation
   follows the dimensionally consistent weak form Eq.27.
2. Reevaluate the UNWEIGHTED residual and strain. Compute the full-tensor
   Eq.16 driving force with Eq.34 block averages:
   `Phi=w1*tanh((mean(f)*mean(r)^2*(epsilon:epsilon)-a1)/c1)
        +w2*tanh((mean(f)-a2)/c2)`.
   Signed r is averaged over non-overlapping 4x4 blocks before squaring;
   f is averaged over separate non-overlapping 8x8 blocks. This is not a
   moving-window mean-absolute-residual screening. f here is the normalized
   reference image at the fine-grid locations, as used in the paper's formula.
3. With u fixed solve Eq.33 for independent nodal damage:
   `[lc*K + diag(mass*(mu/lc+r^2+Phi))]d = mass*(r^2+Phi)`
   subject to the physical bounds 0<=d<=1. The homogeneous natural boundary
   condition is zero normal damage flux. There is no prescribed d=1 crack,
   no distance-profile ansatz, and no imposed irreversibility.
4. Save the state and check the Section 2.7 criterion
   `||M(d)*delta_u||_2 <= tolerance*||M(d)*u_1||_2`.
   Pixel norms include quadrature mass. Damage KKT convergence is also
   required; a stalled displacement solve is not reported as convergence.

## Discretization and numerical safeguards

Displacement uses a structured Q1 mesh with 16 px spacing; damage uses unit
pixel spacing with two linear triangles per pixel cell. The diffusion stiffness
is exact for these triangles, and reaction/load integration uses vertex-lumped
mass. Coarse FE shape functions evaluate displacement and strain on the damage
grid. Only complete cells enter the damage mesh. These integration/element
choices are implementation details, not an assertion of identical author code.

Both images are normalized together to [0,1]. The default inverts them because
these experimental cracks are dark, whereas the paper assumes bright
crack/background pixels. The original uninverted image is used for plotting.

The signed tanh force is kept intact. It can make the damage Hessian indefinite;
silently replacing Phi or r^2+Phi with a positive part would change the model.
A bound-constrained active-set Newton method with red/black coordinate-minimization sweeps and projected-gradient fallback
instead solves the discrete quadratic and records its projected KKT residual.
Bounds, solver tolerances, damping and line searches are explicit numerical
choices because the paper does not specify their implementation.

For these images, a direct one-step zero-start GN run began fitting speckle
mismatches with diffuse damage. To improve the registration capture range,
the first fixed-damage displacement subproblem defaults to Gaussian scales
[4 2 1 0] px, always ending on the original normalized images, with up to
12 inner GN iterations per scale. Subsequent outer steps use only original
images. This coarse-to-fine inner solve is an extension beyond the illustrative
one-step algorithm, not a different damage law. The final original-image
weighted objective must decrease. To run the direct version:

```matlab
p.displacementSmoothing = 0;
p.innerMaxIterations = 1;
```

Set these immediately after `p=paperPFParameters` in the script. Initialization
remains zero in both versions. Defaults also include a 2 px GN step cap and
small diagonal damping; every accepted step records its objective change.

## Parameters and calibration

Defaults: mu=0.1, lc=5 px, displacement spacing 16 px, damage spacing 1 px,
and outer tolerance 0.05. Eq.16 parameters start from the paper's experimental
Table 3:

| Parameter | w1 | a1 | c1 | w2 | a2 | c2 |
|---|---:|---:|---:|---:|---:|---:|
| Value | 0.05 | 3e-5 | 1e-5 | 0.04 | 0.85 | 0.05 |

These were calibrated for the authors' specimen, not this one. Change them
in `paperPFParameters` only with a documented baseline/calibration. Numerical
convergence is not proof of accurate crack geometry. Indicator maxima, damage
support, residual RMSE, and both subsolver diagnostics are saved for inspection.
No ground-truth accuracy percentage is claimed.

## Saved damage and crack overlays

Outputs go to `PaperPhaseField/Run_<timestamp>`:
- `Model.mat`: both mesh operators, normalization and plot bounds.
- `Iteration_00.mat/png`: zero-initialized state.
- `Iteration_XX.mat/png`: every completed displacement/damage iteration.
- `Damage.mat`: final state, also copied to `PaperPhaseField/Damage.mat`.
- `Convergence.mat`: actual status and full iteration history.
- `PaperPhaseField/LatestRun.mat`: pointer to the latest completed run.

The PNGs compare the original deformed crack region with d on a fixed [0,1]
scale and show d=0.5 contours (cyan) and the skeleton (red) over the original
grayscale image. The existing extracted crack may set ONLY the plot zoom.
The full ROI damage is saved regardless of zoom.

The damage>=0.5 contour follows the paper's experimental comparison. Skeleton
and optional principal path are visualization postprocessing, not an extra
constraint in PF-DIC. All branches remain in `crackSkeleton`; the optional
`crackPathDeformedXY` is the longest retained connected path, with no gap
bridges. No path is fabricated if the threshold produces no crack.

Snapshots save Damage, raw normalized residual, driving force, Eulerian nodal
and pixel displacement, both crack representations, parameters and history.
Outside-ROI damage stays zero. These outputs do not overwrite the older
Correli initialization, damage or level-set results.

## Tests

Run the following:

```matlab
results = runtests('tests');
assertSuccess(results);
```

 Tests cover partition
of unity, exact affine coarse-mesh interpolation, diffusion nullspace and mass,
warp sign/Jacobian finite differences, signed non-overlapping averages, Eq.16,
constant-source damage, damage bounds/diffusion, and genuine weighted
nonlinear displacement registration.

## Results on this repository's images

All 43 tests pass. The Table 3 baseline completed two outer iterations with
maximum damage 0.4848, below the unchanged d=0.5 contour threshold; its
crack overlay is therefore empty. This is a saved baseline, not a successful
crack measurement.

Three fixed-displacement parameter trials were then compared visually. The
most restrictive trial was rerun as a complete coupled solve with these
explicit provisional overrides:

```matlab
paperParameterOverrides = struct( ...
    'w1',0.5, ...
    'a1',3e-6, ...
    'c1',3e-7);
runPaperPhaseFieldDIC
```

Everything else stays at the default values. The coupled run converged in
four outer iterations; its weighted increment ratio was 0.04719, maximum
damage 0.7842, and 1646 ROI pixels reached d>=0.5. The damage subproblem
projected KKT norm was approximately 4.1e-15. The damage follows much of the
visible crack, but the 0.5 contour remains fragmented; the centerline is not
a complete validated crack path. No postprocessing bridges were inserted.

The change from the paper's Table 3 parameters is a provisional specimen
calibration, not an author parameter set or independent validation. It is
recorded in every MAT file. To return to the Table 3 baseline:

```matlab
clear paperParameterOverrides
runPaperPhaseFieldDIC
```

The zero-start direct one-step diagnostic run and an earlier nonconverged
damage-solver trial are not used as final results. Fixed-displacement
calibration figures are explicitly labeled as such rather than as staggered
iterations. Neither parameter calibration nor numerical convergence alone
establishes physical crack-location accuracy.
