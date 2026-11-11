function [deformedXY,displacementXY] = mapCrackPathToDeformed(referenceXY,Mesh,a)
%MAPCRACKPATHTODEFORMED Map reference crack nodes by x = X + u(X).
% The residual f(X)-g(X+u(X)) and level sets are on the reference grid.
% Interpolate with the actual Correli T3 elements, not a new triangulation:
% this preserves the mesh's holes and element-wise linear displacement.
% The original reference path is never modified.
validateattributes(referenceXY,{'numeric'}, ...
    {'real','finite','2d','ncols',2,'nonempty'});
validateattributes(Mesh.pos,{'numeric'},{'real','finite','2d','nonempty'});
validateattributes(a,{'numeric'},{'real','finite','2d','nonempty'});
if size(Mesh.pos,2)<2 || size(a,1)~=size(Mesh.pos,1) || size(a,2)<2
    error('mapCrackPathToDeformed:Displacement', ...
        'Displacements must contain ux,uy for every reference-mesh node.');
end
if size(Mesh.con,2)~=4 || isempty(Mesh.con) || any(Mesh.con(:,1)~=0)
    error('mapCrackPathToDeformed:Mesh', ...
        'Expected Correli T3 elements: Mesh.con = [0 node1 node2 node3].');
end
connectivity = double(Mesh.con(:,2:4));
validateattributes(connectivity,{'numeric'}, ...
    {'finite','integer','positive','<=',size(Mesh.pos,1)});
meshTriangles = triangulation(connectivity,double(Mesh.pos(:,1:2)));
[element,barycentric] = pointLocation(meshTriangles,double(referenceXY));
if any(isnan(element))
    error('mapCrackPathToDeformed:OutsideMesh', ...
        '%d crack nodes are outside the displacement mesh; no extrapolation is performed.', ...
        nnz(isnan(element)));
end
vertices = connectivity(element,:);
displacementXY = zeros(size(referenceXY));
for component = 1:2
    values = double(a(:,component));
    displacementXY(:,component) = sum(barycentric.*values(vertices),2);
end
deformedXY = double(referenceXY)+displacementXY;
end
