# Geometry-constrained phase-field DIC

## Run order and saved data

Configure your existing Correli installation so `correlate` and `correli_pathdef` are available. Run from the project root:

```matlab
example2Dscreening
extractResidualCrackPath
runPhaseFieldDIC
```

- `example2Dscreening` now saves `InitialValues/DisplacementInit.mat`: nodal displacements `a`, the translated `Mesh`, original-image residual `R`, `measurementMask=isfinite(R)`, `imageSize`, and the fixed reference gray-level range `grayScale`. It also saves `Residual/R.mat`. The correlator option is corrected to `conv_limit`.
- `extractResidualCrackPath` saves `CrackPaths/CrackPathInit.mat`: ordered `crackPathXY` ([column,row] pixel coordinates), raster mask, image size, and graph diagnostics. Empty paths are not saved.
- `runPhaseFieldDIC` saves `Damage/Damage.mat`: image-sized `Damage`, the fixed measurement mask, initial/current `levelSets`, parameters, energy history and acceptance diagnostics.
- `Damage/PhaseFieldResult.mat` saves accepted `a`, original-image `R`, mesh, parameters and energy history.

The scripts create the folders as needed. Running them replaces their named result files. MAT results are generated locally; only source and empty-folder placeholders are versioned.

## Three level sets

`LS=crackLevelSets2D(pathXY,imageSize)` extends the crack-surface/front construction of Rethore, Hild and Roux, *Extended digital image correlation with crack shape optimization*, DOI 10.1002/nme.2070, Section 4.1, Eqs. (20)-(22), to two tips:

- `LS.phi`: signed Euclidean distance to the oriented polyline extended beyond both ends by tangent rays, zero on the crack surface. Positive is left of the directed curve in image coordinates (y increases downwards).
- `LS.psiStart=-(X-Pstart).tStart` and `LS.psiEnd=(X-Pend).tEnd`: signed distances to lines normal to the endpoint tangents. Their zero contours intersect phi=0 at the respective tips. Unit gradients are orthogonal to the crack-surface gradient locally at each tip.
- `LS.distance`: exact unsigned distance to the *finite* polyline. This avoids treating the tangent extensions as damaged material beyond the tips.
- `LS.pathXY` stores the ordered geometry.

These piecewise distance fields have the usual closest-point/cut-locus nondifferentiabilities. Front planes are local tip coordinates, not a global chart for arbitrarily folded cracks. For a folded curve a front plane can intersect other parts of the crack; the finite-polyline distance remains the authoritative damage geometry. No X-FEM enrichment is implemented in this change.

## Damage and energy

Zhu and Lu, *Phase-field digital image correlation for integrated displacement and damage measurements*, Eq. (8), is used with the extra damage-driving/mechanical term Phi=0:

```text
E = 1/2 sum_ROI (1-d)^2 (R/grayScale)^2
  + 1/2 sum_ROI (mu/lc) d^2
  + lc/2 sum_internal_pixel_bonds (d_j-d_i)^2.
```

Pixel spacing is one in both directions; lengths are in pixels. Only horizontal/vertical bonds with both endpoints in the fixed measurement domain contribute. This gives a natural no-flux ROI boundary, including holes. Missing residuals inside the fixed domain cause an error rather than silently shrinking the domain. Residuals are normalized by a fixed measured reference-image range; this scale and the parameters must be calibrated for the experiment.

`createDamageField` uses the geometry-constrained profile

```text
d = exp[-sqrt(mu)/lc * (sqrt(distance^2 + epsilon^2) - epsilon)].
```

As epsilon tends to zero, this is the exact homogeneous straight-crack profile for `(mu/lc)d-lc*Delta(d)=0`, with d=1 on the crack and exponential decay away from it. Setting mu=1 recovers the profile associated with the paper's Eq. (5). For curved finite cracks this is a distance-based *ansatz*, not an exact PDE solution. Default epsilon=0.1 px smooths its derivative at zero distance; the default lc=5 and mu=0.1 follow the paper's example.

This implements a constrained crack-geometry optimization. It does **not** solve the paper's independent residual-driven damage PDE (Eq. (13)), introduce new cracks, enforce damage irreversibility, or use the additional strain-based driving force. Damage may increase or decrease as the crack is repositioned.

Only ROI pixels are updated. An optional baseDamage array is preserved exactly outside the ROI; the default is zero there. Outside-ROI reference/deformed image pixels remain unchanged. Damage is a dimensionless field, not an image containing gray-level intensities.

Example:
```matlab
parameters=phaseFieldParameters(struct('lc',5,'mu',0.1,'grayScale',grayScale));
LS=crackLevelSets2D(crackPathXY,size(R));
Damage=createDamageField(LS,measurementMask,parameters,[], ...
    fullfile('Damage','Damage.mat'));
```

## Level-set updates

`updateCrackLevelSets` analytically differentiates the discrete energy with respect to d, then pulls that derivative back to polyline nodes through the nearest-segment projection and damage profile. Interior nodes move normally; both endpoint nodes can move in two dimensions. The default maximum movement is 0.25 px. Backtracking rejects non-decreasing energy, collapsed segments, self-intersections, or a path leaving the ROI. Both endpoints and all three distance fields are reconstructed after every accepted move.

This geometric update/reinitialization is an explicit adaptation of the level-set representation; it does **not** implement Rethore's Godunov time integrator or residual-maximum velocity rule. Reconstructing signed distances keeps the geometry and the three fields synchronized. Acceptance guarantees energy decrease only for each accepted discrete step, not a global optimum or convergence after the default three outer iterations.

## Existing correlate integration: approximation and acceptance

`applyDamageToImages` multiplies both images by **1-d**, the square root of the degradation function M=(1-d)^2, inside the ROI. Multiplying by M itself would incorrectly square the intended energy weight.

With a nonzero displacement, fixed-grid image masking is generally **not exactly** the PF-DIC functional:
```text
masked residual: (1-d(X))*f(X) - (1-d(X+u))*g(X+u)
desired residual: (1-d(X)) * [f(X)-g(X+u)].
```
The masked images also change the gradients used by correlate. Thus this is the requested image-only surrogate, not an exact weighted DIC solver. A fully variational displacement solve would need material-coordinate weighting inside the correlator and its linearization.

`correlatePhaseField` calls your existing `correlate` with those images, using the saved nodal displacement as initialization. Mechanical regularization and brightness/contrast correction are disabled. It then uses the returned Correli kernel handle to evaluate the proposal on the **original images**. Only a non-increasing actual PF energy on the same ROI is accepted. Invalid ROI coverage or increased energy causes rejection and preserves the previous displacements/residual. Consequently, it is valid for every masked-image displacement proposal to be rejected while geometry updates proceed.

## Validation

Tests cover straight/curved subpixel distance geometry, tip orthogonality, damage bounds and outside-domain preservation, MAT serialization, finite-difference checks of both damage and shape energy derivatives, energy-decreasing geometry motion, endpoint-field reconstruction, unchanged outside image pixels, and invalid domains/paths.

The full workflow was also run with the installed Correli on this repository's images and mesh. It generated all requested MAT files. Three accepted geometry steps reduced total energy from 248.965917 to 247.641660. All three image-masked displacement proposals were rejected by the actual-energy check, so the regular-DIC displacement remained unchanged.

```matlab
results=runtests('tests');
assertSuccess(results);
```
