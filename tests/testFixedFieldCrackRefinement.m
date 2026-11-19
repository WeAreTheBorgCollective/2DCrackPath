function testFixedFieldCrackRefinement
[x,y] = meshgrid(1:60);
field = 20*exp(-((y-30).^2)/30);
Mesh.pos = [2 2 0;59 2 0;2 59 0;59 59 0];
Mesh.con = [0 1 2 3;0 2 4 3];
path = [10 25;50 25];
original = field;
r = refineCrackPathFixedField(path,field,Mesh,struct('maxIterations',8));
assert(isequal(r.fixedField,original));
assert(numel(r.scoreHistory)>1 && all(diff(r.scoreHistory)>0));
assert(median(r.pathXY(:,2))>25);
assert(crackPathInDomain(r.pathXY,isfinite(field)));
assert(isequal([r.graph.Nodes.X r.graph.Nodes.Y],r.pathXY));
assert(isequal(r.levelSets.pathXY,r.pathXY));
r = refineCrackPathFixedField(path,ones(size(field)),Mesh,struct('maxIterations',2));
assert(isequal(r.pathXY,path) && numel(r.scoreHistory)==1);
fprintf('Fixed field, monotone ascent, contour movement, graph and stationarity checks passed.\n');
end
