function testCrackTipSnapping
Mesh.pos = [5 5 0;45 5 0;5 45 0;45 45 0];
Mesh.con = [0 1 2 3;0 2 4 3];
% Along y=15 the first triangle spans x=5..35, so its directional length is 30.
path = [3 15;33 15];
[q,info] = snapCrackTipsToElements(path,Mesh);
assert(norm(q(end,:)-[35 15])<1e-10);
assert(info(2).snapped && abs(info(2).remainingDistance-2)<1e-10);
assert(abs(info(2).directionalLength-30)<1e-10);
% Exactly one quarter away does not satisfy the strict threshold.
[q,info] = snapCrackTipsToElements([3 15;27.5 15],Mesh);
assert(~info(2).snapped && q(end,1)==27.5);
[q,info] = snapCrackTipsToElements([3 15;20 15],Mesh);
assert(~info(2).snapped && q(end,1)==20);
% Border tips do not cascade into the adjacent element on repeated calls.
[q,~] = snapCrackTipsToElements(path,Mesh);
[again,info] = snapCrackTipsToElements(q,Mesh);
assert(isequal(q,again) && ~any([info.snapped]));
[reverse,~] = snapCrackTipsToElements(flipud(path),Mesh);
assert(max(abs(reverse(:)-reshape(flipud(q),[],1)))<1e-10);
% A crack wholly inside one element does not count as an almost completed cut.
[q,info] = snapCrackTipsToElements([30 15;33 15],Mesh);
assert(~any([info.snapped]) && isequal(q,[30 15;33 15]));
% A forbidden target is left unchanged.
mask = true(50);
mask(:,35:end) = false;
[q,info] = snapCrackTipsToElements(path,Mesh,struct,mask);
assert(~info(2).snapped && isequal(q,path));
% MEX: tip interior -> partial cut; exit-edge tip -> full cut of ONLY triangle 1.
partial = prepareXDIC(Mesh,[50 50],path);
edge = prepareXDIC(Mesh,[50 50],[3 15;35 15]);
assert(isempty(partial.fullyCutElements));
assert(isequal(edge.fullyCutElements,1));
assert(~ismember(2,edge.cutElements));
assert(~isempty(edge.enrichedNodes));
% Crack running along an element edge does not cross either interior.
boundary = prepareXDIC(Mesh,[50 50],[5 45;45 5]);
assert(isempty(boundary.cutElements));
fprintf('Tip snapping, strict quarter threshold, repeat stability and MEX border classification passed.\n');
end
