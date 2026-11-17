function testXDICROIBoundary
% A fractional crack endpoint belongs to the ROI but has an outside neighbor.
[x,y] = meshgrid(1:60);
mask = x>=10 & x<=50 & y>=5 & y<=55;
R = double(y);
R(~mask) = NaN;
q = [9.8 20.25;10.2 20.25;8.8 20.25];
v = sampleROIImage(R,mask,q);
assert(max(abs(v(1:2)-20.25))<1e-12);
assert(isnan(v(3)));
[gx,gy] = gradientInROI(R,mask);
assert(all(abs(gx(mask))<1e-12));
assert(all(abs(gy(mask)-1)<1e-12));
% Previously pathScore returned -Inf at the first endpoint and froze phi.
LS = crackLevelSets2D([9.8 20;40 20],size(R));
[new,info] = updateXDICLevelSets(LS,R,mask, ...
    struct('screeningWindowSize',[1 1],'gradientGain',0.5,'smoothLength',0));
assert(info.accepted && isfinite(info.meanScreeningBefore));
assert(info.meanScreeningAfter>info.meanScreeningBefore);
assert(median(new.pathXY(:,2))>20.2);
assert(crackPathInDomain(new.pathXY,mask));
% A constant field stays stationary, including at the ROI border.
R(mask) = 7;
[unchanged,info] = updateXDICLevelSets(LS,R,mask,struct);
assert(~info.accepted && isequal(unchanged.pathXY,LS.pathXY));
assert(isfinite(info.meanScreeningBefore));
Mesh.pos = [5 5 0;45 5 0;5 45 0;45 45 0];
Mesh.con = [0 1 2 3;0 2 4 3];
[constrained,detail] = constrainCrackToMesh([4.8 15;20 20;40 15],Mesh);
tri = triangulation(Mesh.con(:,2:4),Mesh.pos(:,1:2));
assert(detail.valid && detail.projectedPoints==1);
assert(all(isfinite(pointLocation(tri,constrained))));
assert(norm(constrained(1,:)-[5 15])<1e-6);
[~,detail] = constrainCrackToMesh([2 15;20 20],Mesh);
assert(~detail.valid);
fprintf('ROI-aware interpolation, one-sided gradients and boundary contour motion passed.\n');
end
