# Snap nearly completed triangle cuts

XDIC evolution now snaps tips before the initial solve and before re-solving
shape and face-growth proposals. options.tipSnap.enabled defaults to true
and fraction defaults to 0.25.

For a tip strictly inside a triangle, the last nonzero crack segment defines
its outward direction. Intersect the line through that tip with the triangle:
L is the full entry-to-exit chord length in that direction, and d is the
remaining forward distance to the exit edge. Snap only when d < 0.25*L.
This defines "element length" directionally, rather than as its longest edge.
The existing terminal segment is extended collinearly to the exact exit.
Both endpoints are treated independently, and the crack must actually enter
the triangle (a crack wholly contained in one element is not a completed cut).
An endpoint already on an edge/vertex is unchanged, preventing repeated
snapping through neighboring elements.

The MEX preparation kernel now classifies only strictly interior tips as
tip-containing. A fully traversed triangle with an endpoint on its edge is
eligible for Heaviside enrichment. A segment touching only an edge or vertex
does not mark the neighboring triangle as cut. No singular tip enrichments
are added. Rebuild with buildXDICKernels after updating the source.

Within evolution, domain and non-crossing checks still apply. Changes to
shape/growth proposals retain the existing re-solve and residual acceptance
test; the initial snap establishes the starting geometry. Standard and
enriched warm starts are transferred by node ID when the basis changes.
All level sets, paths and saved bases use the same snapped geometry.

The 10 px limit still caps the face-intersection proposal. This new element
completion step is additional and may move the tip farther than that cap.
Snapshot info.tipSnap records element ID, chord length, remaining distance,
extension, target and rejection reason; result.options stores the settings.
Outside-mesh, already-on-border and invalid/domain-crossing snaps are skipped.

Validation: snapping/strict-quarter threshold, reversed ordering, repeated
border calls, wholly-contained cracks and masked targets passed. Rebuilt MEX
tests confirm the traversed triangle becomes fully cut, the touching neighbor
does not, and edge-aligned cracks do not cut triangle interiors. Existing
basis, Hessian/RHS, jump recovery, face-growth and level-set tests passed.
An isolated evolution run verified initial snapping, enriched warm-start
transfer, level sets and saved basis agreement. The currently saved original
initial crack does not satisfy the quarter-length threshold at either end,
so it is correctly unchanged by the rule alone.
