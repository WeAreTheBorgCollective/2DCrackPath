# Compare initial and final phase-field DIC

Run `comparePhaseFieldDIC` after `runPhaseFieldDIC`, with Correli configured on the MATLAB path.

The script reads the initial displacement from `InitialValues/DisplacementInit.mat` and the **final accepted** displacement from `Damage/PhaseFieldResult.mat`. It recomputes both residuals with the native Correli cubic kernel on the original reference and deformed images. No damage weighting is applied and no correlation optimization is rerun.

The kernel samples the original deformed image at X+u(X), giving R=f(X)-g(X+u(X)). The displayed warped image is reconstructed as f-R, retaining the kernel's interpolation and mesh mapping. Floating-point residual storage introduces only roundoff; an independent zero-reference kernel evaluation agreed within 0.001954 GL.

Figures use a common valid ROI and shared initial/final color limits. Coverage counts are reported separately; losses are warned about. Displacement plots interpolate nodal x/y components over the actual triangular FE mesh, not over a convex-hull approximation.

Outputs:
- `InitialValues/DisplacementInitPF.mat`: final nodal displacement under both `a` and `DisplacementInitPF`, mesh, recomputed unmasked `R`, measurement/comparison masks, image size, gray-level scale, and comparison metrics.
- `Comparisons/ResidualComparison.png`: reference, initially/finally warped deformed image, signed residuals, and their difference.
- `Comparisons/DisplacementComparison.png`: initial/final/difference maps for x and y displacement components.
- `Comparisons/PhaseFieldComparison.mat`: warped images, residuals, common mask, metrics and nodal changes.

All outputs are generated locally when the script runs. Existing named outputs are replaced.

For the validated saved run, the two displacement arrays were identical because all masked-image displacement proposals had been rejected by the phase-field energy check. Both unmasked residual RMSE values are 122.639 GL, MAE values are 82.522 GL, and the common ROI contains 407,035 pixels. The phase-field energy decrease in that run came from geometry/damage changes, not improved displacement or unmasked residual. The comparison does not imply a displacement improvement.
