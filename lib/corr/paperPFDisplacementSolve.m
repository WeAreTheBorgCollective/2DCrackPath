function [updated,increment,r,info] = paperPFDisplacementSolve(model,f,g,U,d,p,outerIteration)
%PAPERPFDISPLACEMENTSOLVE Fixed-damage nonlinear weighted registration.
% Coarse-to-fine Gaussian images on the FIRST outer step avoid speckle-scale
% local minima with zero initialization. This is a documented numerical
% extension, not a change to the damage functional. Last scale is always zero.
% Set displacementSmoothing=0 and innerMaxIterations=1 for the direct version.
start = U;
if outerIteration==1
    scales = p.displacementSmoothing;
else
    scales = 0;
end
if isempty(scales) || scales(end)~=0
    error('paperPF:Scales','Displacement smoothing must finish at the original images (0).');
end
weights = model.mass.*(1-d).^2;
startR = g(model.ids)-paperPFWarp(f,model.x,model.y,model.N*start);
oldCost = 0.5*sum(weights.*startR.^2);
count = 0;
lastInfo = struct('gradientNorm',Inf);
for sigma = scales
    if sigma>0
        fs = imgaussfilt(f,sigma,'Padding','replicate');
        gs = imgaussfilt(g,sigma,'Padding','replicate');
    else
        fs = f;
        gs = g;
    end
    for inner = 1:p.innerMaxIterations
        [trial,delta,~,lastInfo] = paperPFDisplacementStep(model,fs,gs,U,d,p);
        count = count+1;
        if ~lastInfo.accepted
            break
        end
        U = trial;
        if norm(delta,'fro')<=p.innerTolerance*max(1,norm(U,'fro'))
            break
        end
    end
end
r = g(model.ids)-paperPFWarp(f,model.x,model.y,model.N*U);
newCost = 0.5*sum(weights.*r.^2);
if newCost>oldCost || any(~isfinite(r))
    [U,~,r,lastInfo] = paperPFDisplacementStep(model,f,g,start,d,p);
    newCost = 0.5*sum(weights.*r.^2);
end
updated = U;
increment = U-start;
info = struct('accepted',norm(increment,'fro')>0,'innerIterations',count, ...
    'smoothingScales',scales,'weightedCostBefore',oldCost, ...
    'weightedCostAfter',newCost,'gradientNorm',lastInfo.gradientNorm);
end
