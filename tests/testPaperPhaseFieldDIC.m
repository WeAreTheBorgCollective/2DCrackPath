function tests = testPaperPhaseFieldDIC
tests = functiontests(localfunctions);
end

function setupOnce(~)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setupPaths;
end

function testTwoMeshOperators(tc)
mask = true(13,17);
model = paperPFModel(mask,4);
verifyEqual(tc,full(sum(model.N,2)),ones(nnz(mask),1),'AbsTol',1e-12);
verifyEqual(tc,model.N*model.nodes,[model.x model.y],'AbsTol',1e-12);
verifyEqual(tc,model.Bx*model.nodes(:,1),ones(nnz(mask),1),'AbsTol',1e-12);
verifyEqual(tc,model.By*model.nodes(:,2),ones(nnz(mask),1),'AbsTol',1e-12);
verifyEqual(tc,model.K*ones(nnz(mask),1),zeros(nnz(mask),1),'AbsTol',1e-12);
verifyEqual(tc,sum(model.mass),12*16,'AbsTol',1e-12);
verifyEqual(tc,model.K,model.K','AbsTol',1e-12);
end

function testWarpSignAndDerivative(tc)
[x,y] = meshgrid(1:20,1:20);
f = 2*x+3*y+0.01*x.*y;
points = [7.3 8.2;11.2 13.3];
u = [0.2 -0.1;0.1 0.2];
[warped,gx,gy,valid] = paperPFWarp(f,points(:,1),points(:,2),u);
q = points-u;
verifyTrue(tc,all(valid));
verifyEqual(tc,warped,2*q(:,1)+3*q(:,2)+0.01*q(:,1).*q(:,2),'AbsTol',1e-12);
for component = 1:2
    plus = u;
    minus = u;
    h = 1e-6;
    plus(:,component) = plus(:,component)+h;
    minus(:,component) = minus(:,component)-h;
    rp = -paperPFWarp(f,points(:,1),points(:,2),plus);
    rm = -paperPFWarp(f,points(:,1),points(:,2),minus);
    gradients = [gx gy];
    verifyEqual(tc,(rp-rm)/(2*h),gradients(:,component),'AbsTol',1e-7);
end
end

function testSignedNonoverlappingAverages(tc)
mask = true(4,8);
values = repmat([1 -1 1 -1 2 2 2 2],4,1);
averaged = paperPFBlockMean(values(:),mask,4);
expected = [zeros(4,4) 2*ones(4,4)];
verifyEqual(tc,averaged,expected(:));
end

function testIndependentDamagePDE(tc)
model = paperPFModel(true(15,15),5);
p = paperPFParameters;
H = 0.2*ones(nnz(model.mask),1);
[d,info] = paperPFDamageSolve(model.K,model.mass,H,zeros(size(H)),p,model.colors);
verifyTrue(tc,info.converged);
verifyEqual(tc,d,H./(p.mu/p.lc+H),'AbsTol',1e-7);
H(:) = -0.05;
[d,info] = paperPFDamageSolve(model.K,model.mass,H,zeros(size(H)),p,model.colors);
verifyTrue(tc,info.converged);
verifyEqual(tc,d,zeros(size(H)));
end

function testDamageDiffusionAndBounds(tc)
model = paperPFModel(true(21,21),5);
p = paperPFParameters;
H = zeros(nnz(model.mask),1);
H(221) = 3;
[d,info] = paperPFDamageSolve(model.K,model.mass,H,zeros(size(H)),p,model.colors);
verifyTrue(tc,info.converged);
verifyGreaterThan(tc,d(220),0);
verifyGreaterThan(tc,d(221),d(220));
verifyGreaterThanOrEqual(tc,min(d),0);
verifyLessThanOrEqual(tc,max(d),1);
end

function testTrueWeightedDisplacementStep(tc)
[x,y] = meshgrid(1:50,1:50);
f = 0.5+0.2*sin(x/3)+0.15*cos(y/4);
g = reshape(paperPFWarp(f,x(:),y(:),repmat([0.2 -0.15],numel(x),1)),size(f));
mask = false(size(f));
mask(8:42,8:42) = true;
model = paperPFModel(mask,8);
p = paperPFParameters;
U = zeros(size(model.nodes));
[updated,~,~,info] = paperPFDisplacementStep(model,f,g,U,zeros(nnz(mask),1),p);
verifyTrue(tc,info.accepted);
verifyLessThan(tc,info.weightedCostAfter,info.weightedCostBefore*0.1);
verifyEqual(tc,median(model.N*updated),[0.2 -0.15],'AbsTol',0.02);
end

function testDrivingForceEq16(tc)
model = paperPFModel(true(9,9),4);
U = [0.1*model.nodes(:,1) 0.2*model.nodes(:,2)];
r = 0.3*ones(nnz(model.mask),1);
p = paperPFParameters;
f = 0.8*ones(9);
[Phi,details] = paperPFDrivingForce(model,U,r,f,p);
expectedIndicator = 0.8*0.3^2*(0.1^2+0.2^2);
verifyEqual(tc,details.indicator,expectedIndicator*ones(size(r)),'AbsTol',1e-12);
expectedPhi = p.w1*tanh((expectedIndicator-p.a1)/p.c1)+p.w2*tanh((0.8-p.a2)/p.c2);
verifyEqual(tc,Phi,expectedPhi*ones(size(r)),'AbsTol',1e-12);
end

function testNonlinearWeightedSolve(tc)
[x,y] = meshgrid(1:50,1:50);
reference = 0.5+0.2*sin(x/3)+0.15*cos(y/4);
deformed = reshape(paperPFWarp(reference,x(:),y(:), ...
    repmat([0.2 -0.15],numel(x),1)),size(reference));
mask = false(size(reference));
mask(8:42,8:42) = true;
model = paperPFModel(mask,8);
p = paperPFParameters;
p.displacementSmoothing = 0;
d = 0.3*ones(nnz(mask),1);
[U,delta,r,info] = paperPFDisplacementSolve( ...
    model,reference,deformed,zeros(size(model.nodes)),d,p,1);
verifyTrue(tc,info.accepted);
verifyEqual(tc,U,delta);
verifyLessThan(tc,info.weightedCostAfter,info.weightedCostBefore*1e-5);
verifyEqual(tc,median(model.N*U),[0.2 -0.15],'AbsTol',1e-3);
verifyEqual(tc,info.weightedCostAfter, ...
    0.5*sum(model.mass.*(1-d).^2.*r.^2),'AbsTol',1e-12);
end

function testDamageCanDecreaseAfterNegativeForce(tc)
model = paperPFModel(true(15,15),5);
p = paperPFParameters;
H = -0.5*ones(nnz(model.mask),1);
[d,info] = paperPFDamageSolve(model.K,model.mass,H, ...
    0.9*ones(size(H)),p,model.colors);
verifyTrue(tc,info.converged);
verifyEqual(tc,d,zeros(size(H)),'AbsTol',1e-6);
end
