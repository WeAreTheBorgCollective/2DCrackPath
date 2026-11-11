function testXDICLevelSets
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setupPaths;
[~,y] = meshgrid(1:80);
LS = crackLevelSets2D([(20:60)' 40*ones(41,1)],size(y));
for direction = [-1 1]
    R = 100+direction*y;
    options = struct('screeningWindowSize',[1 1],'gradientGain',.5);
    [updated,info] = updateXDICLevelSets(LS,R,true(size(R)),options);
    assert(info.accepted);
    assert(abs(median(updated.pathXY(:,2))-40-direction*.5)<.05);
    assert(abs(interp2(updated.psiStart,updated.pathXY(1,1),updated.pathXY(1,2)))<1e-10);
    assert(abs(interp2(updated.psiEnd,updated.pathXY(end,1),updated.pathXY(end,2)))<1e-10);
    assert(max(abs(info.velocity))<=1);
end
[updated,info] = updateXDICLevelSets(LS,ones(80),true(80));
assert(~info.accepted);
assert(isequal(updated,LS));
mask = true(80);
mask(43:44,:) = false;
R = 100+y;
R(~mask) = NaN;
[updated,~] = updateXDICLevelSets(LS,R,mask,struct('screeningWindowSize',[1 1]));
assert(crackPathInDomain(updated.pathXY,mask));
fprintf('PASS: gradient directions, front intersections, speed bound, flat residual and domain mask.\n');
end
