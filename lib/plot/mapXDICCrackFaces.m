function [minusFace,plusFace] = mapXDICCrackFaces(path,Mesh,solution,prep)
%MAPXDICCRACKFACES Map both one-sided displacement traces to the deformed image.
% Averaging these traces is a display centerline, not a single crack displacement.
tri = triangulation(double(Mesh.con(:,2:4)),double(Mesh.pos(:,1:2)));
[element,bary] = pointLocation(tri,path);
if any(isnan(element))
    error('xdic:Mapping','Crack point outside the actual mesh.');
end
nodes = Mesh.con(element,2:4);
amplitudes = zeros(size(solution.standard));
amplitudes(solution.enrichedNodes,:) = solution.enriched;
minusFace = path;
plusFace = path;
for component = 1:2
    standard = solution.standard(:,component);
    enriched = amplitudes(:,component);
    minusValue = standard-enriched.*prep.nodeHeaviside;
    plusValue = standard+enriched.*(1-prep.nodeHeaviside);
    minusFace(:,component) = path(:,component)+sum(bary.*minusValue(nodes),2);
    plusFace(:,component) = path(:,component)+sum(bary.*plusValue(nodes),2);
end
end
