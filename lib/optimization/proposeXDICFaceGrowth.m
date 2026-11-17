function [path,info] = proposeXDICFaceGrowth(LS,Mesh,solution,prep,options)
% Build physical one-sided traces and local mean-trace deformation gradients.
[minusFace,plusFace] = mapXDICCrackFaces(LS.pathXY,Mesh,solution,prep);
tri = triangulation(double(Mesh.con(:,2:4)),double(Mesh.pos(:,1:2)));
elements = pointLocation(tri,LS.pathXY([1 end],:));
assert(all(isfinite(elements)),'Crack tips must be inside the mesh.');
amplitude = zeros(size(solution.standard));
amplitude(solution.enrichedNodes,:) = solution.enriched;
meanDisplacement = solution.standard+amplitude.*(0.5-prep.nodeHeaviside);
tipF = zeros(2,2,2);
for tip = 1:2
    nodes = Mesh.con(elements(tip),2:4);
    xy = Mesh.pos(nodes,1:2);
    u = meanDisplacement(nodes,:);
    gradient = (xy(2:3,:)-xy(1,:))\(u(2:3,:)-u(1,:));
    tipF(:,:,tip) = eye(2)+gradient';
end
[path,info] = growCrackFromFaceTraces(LS.pathXY,minusFace,plusFace,tipF,options);
end
