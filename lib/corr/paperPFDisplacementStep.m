function [updated,increment,r,info] = paperPFDisplacementStep(model,f,g,U,d,p)
%PAPERPFDISPLACEMENTSTEP Eq.21 weighted Gauss-Newton, damage held fixed.
% Weights live at fixed deformed-grid quadrature points. Neither image is masked.
u = model.N*U;
[warped,gx,gy,valid] = paperPFWarp(f,model.x,model.y,u);
if ~all(valid)
    error('paperPF:Coverage','Current displacement samples outside the reference image.');
end
r = g(model.ids)-warped;
M = (1-d).^2;
weights = model.mass.*M;
J = [spdiags(gx,0,numel(r),numel(r))*model.N ...
     spdiags(gy,0,numel(r),numel(r))*model.N];
Hessian = J'*(spdiags(weights,0,numel(r),numel(r))*J);
gradient = J'*(weights.*r);
scaling = max(full(diag(Hessian)),1e-12);
delta = -(Hessian+p.damping*spdiags(scaling,0,numel(gradient),numel(gradient)))\gradient;
delta = reshape(delta,size(U));
maxStep = max(vecnorm(delta,2,2));
if maxStep>p.maxDisplacementStep
    delta = delta*(p.maxDisplacementStep/maxStep);
end
oldCost = 0.5*sum(weights.*r.^2);
accepted = false;
alpha = 0;
updated = U;
increment = zeros(size(U));
for backtrack = 0:20
    trialAlpha = 2^(-backtrack);
    trial = U+trialAlpha*delta;
    [warped,~,~,valid] = paperPFWarp(f,model.x,model.y,model.N*trial);
    trialR = g(model.ids)-warped;
    if all(valid) && 0.5*sum(weights.*trialR.^2)<oldCost
        updated = trial;
        increment = trialAlpha*delta;
        r = trialR;
        alpha = trialAlpha;
        accepted = true;
        break
    end
end
info = struct('accepted',accepted,'alpha',alpha, ...
    'weightedCostBefore',oldCost,'weightedCostAfter',0.5*sum(weights.*r.^2), ...
    'gradientNorm',norm(gradient,Inf));
end
