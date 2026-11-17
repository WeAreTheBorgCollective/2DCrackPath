# Screening cases and image sequences

Run `example2Dscreening` for the default `composite100` case.
Select `screeningCaseName = 'composite050'`, `'composite100'`,
`'composite200'`, `'windmill'`, or `'original'` before running the script.
Edit `screeningCases.m` to add a case with an ordered image list and mesh.

The uploaded filenames are composite00.bmp through composite11.bmp.
The first image remains the reference for all eleven solves. The preceding
total displacement initializes the next Gaussian-smoothed solve (sigma 2);
the unfiltered solve then starts from that coarse result. Displacements are
not summed and the reference mesh is not moved between steps.

MAT meshes use mesh.pos(:,1:2) and mesh.con(:,2:4). The imported z and element
type are discarded. Correli still requires three coordinate columns and an
element-type column, so zeros are inserted for z and the T3 type. Clockwise
triangles are reordered. No coordinate offset is applied to the new meshes.
The original case retains its previous mesh translation.

Each run has a timestamped directory under Results/Screening/<case>/.
Every Step_NNN contains DisplacementInit.mat, R.mat, Images.png and
Residual.png. The mesh overlay is saved once per run as Mesh.png. Figures
retain the previous reference/deformed pair, mesh overlay and signed
grayscale residual displays, with the same windows reused for each step.

After every successful solve, InitialValues/DisplacementInit.mat and
Residual/R.mat are refreshed with that step. Thus they contain step 11
after a complete composite run. All original variables remain available;
image filenames, case information, initial guess and solver statuses are
also saved. A finite result does not by itself prove solver convergence:
inspect status and coarseStatus when assessing a run.

The extraction, standard XDIC and phase-field entry points now read the
active image filenames from DisplacementInit.mat; old files retain the
reference.tif/deformed.tif fallback. The dedicated runXDICCrackMesh2 example
remains tied to its own original dataset. Regenerate crack/damage outputs
after switching cases or timesteps; an old crack file is not automatically
valid for a new pair.

To activate an archived timestep, copy its DisplacementInit.mat into
InitialValues and its R.mat into Residual, then rerun extraction.
For batch use, call runScreeningCase(projectRoot,cases.composite100).
Options ShowPlots=false and PublishLatest=false disable plots or updating
the active legacy files respectively, while preserving per-step MAT results.

Validation (2026-11-17): all three uploaded meshes imported successfully
(506/206/116 nodes for 050/100/200). All eleven composite100 steps completed
with both solver stages reporting convergence. Archived initial guesses
exactly matched preceding solutions; residual files, ROI masks and image
metadata passed regression checks. The original mesh and legacy image
fallback also passed. The validation run used PublishLatest=false to retain
the pre-existing active initialization. Generated results are local, not
committed to the repository.

## Windmill experiment (fifth case)

Select `screeningCaseName = 'windmill'; example2Dscreening`.
The ordered list uses Windmill/windimill_meso_spec_1_001.tif as the fixed
reference, followed by acquisition IDs 050:50:950 and 960 (21 images,
20 displacement solutions). Mesh/Mesh.mat is imported with zero offset.
The spelling windimill matches the uploaded filenames.

The existing solver, Gaussian initialization, previous-solution warm starts
and plots are reused. Results go to Results/Screening/windmill/Run_*/Step_001
through Step_020. Step_020 corresponds to acquisition 960; the image filenames
saved in each MAT identify the acquisition explicitly. Normal runs update
InitialValues/DisplacementInit.mat and Residual/R.mat at each step for the
existing downstream workflows.

Windmill validation (2026-11-17): all 20 timesteps converged in both solver
stages. Previous-solution warm starts, image metadata, residual files and
all per-step plots were verified. Validation used PublishLatest=false;
the previously active initialization was preserved.
