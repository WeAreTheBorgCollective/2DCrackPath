# Residual crack path extraction

Run `extractResidualCrackPath` from the project root after `screenResidualImage` to reuse `rho`. It also runs standalone: if `rho` is absent, it loads `Residual/R.mat` and computes an 11 x 11 mean-absolute-residual map. Clear `rho` to recompute with changed window settings. Requires MATLAB and Image Processing Toolbox (`bwskel`).

## Thresholding
`[BW,T,info] = valleyOtsuThreshold(rho,256)` implements Ng, *Automatic thresholding for defect detection* (2006), Section 2.2, Eq. (8), DOI: 10.1016/j.patrec.2006.03.009:
`(1-p(t)) * (w1(t)*mu1(t)^2 + w2(t)*mu2(t)^2)`.
This is not the same as multiplying between-class variance by the valley weight. Scores are nonnegative. Finite data alone form the histogram; invalid pixels remain background. Values are quantized by rounding to 256 uniformly spaced levels between zero and the finite maximum. Zero is retained because subtracting a data-dependent minimum changes the weighted objective. The reported threshold is the upper boundary of the chosen bin in original GL units; foreground is `rho >= T`. Constant/empty finite data produce no foreground. Equal maxima use the first bin. Histogram bin count is configurable and can affect the result.

## Skeleton and graph
`thinCrackMask(BW)` produces a one-pixel-wide skeleton with MATLAB's topology-preserving `bwskel`, using a false border. This is an explicit implementation substitution for the paper's four-subiteration thinning routine; it is not a reproduction of that routine.

`[xy,mask,G,info] = connectCrackGraph(skeleton,options)` adapts Strini and Schiavi, *Euclidean Graphs as Crack Pattern Descriptors for Automated Crack Analysis in Digital Images* (2022), DOI: 10.3390/s22165942, Sections 3.2, 3.3.1, 3.3.2 and 3.3.5:
- Every skeleton pixel is a graph node, with x=column and y=row. Every 8-neighbor connection has Euclidean weight 1 or sqrt(2).
- Components smaller than the endpoint diameter limit, or exceeding the vertex-count/diameter limit, are removed. Defaults 25 px and 2 vertices/px follow the paper's example, but are not calibrated to this data. Components without two endpoints use all vertices to compute diameter.
- The nearest feasible endpoints in different components are joined by a distal edge, then components and endpoints are recomputed. This conservative, cycle-free merging policy differs from indiscriminately inserting every eligible link. The default maximum gap is 100 px; `Inf` explicitly removes the distance cap. Bridges may not cross false pixels in `ValidMask` (the script uses finite `rho`).
- Weighted shortest paths select the terminal pair with the greatest geodesic distance. Only this principal path is returned, removing side branches and cycles from the path output while retaining them in `G`. With multiple remaining components, the longest principal path is selected and the script warns. With fewer than two endpoints, all vertices are considered as terminals. This fallback can be expensive for large loops.
- This implements the requested main-path subset, not the full paper pipeline: subgraph recovery, secondary-branch recovery and cross-crack correlation are not implemented.

The returned `xy` is an ordered polyline of graph nodes; long segments are inferred gap connections, not observed crack pixels. `mask` rasterizes these segments. `G.Edges.IsBridge` distinguishes inferred edges, plotted in red. `info` reports removed and remaining component counts, excluded nodes, bridge count and Euclidean path length. A single principal path does not mean every retained component has been connected. No path is invented when no component survives filtering.

The default script retains `binaryCrack`, `automaticThreshold`, `crackSkeleton`, `crackGraph`, `crackPathXY`, `crackPathMask`, and diagnostics in the workspace. It does not write or modify input data.

Run tests with `results=runtests('tests'); assertSuccess(results);`. Tests cover the exact Ng objective, quantization/invalid data, thinning connectivity, fragmented-path reconstruction, gap limits, invalid-region barriers, Euclidean distances, branching, loops and empty input.
