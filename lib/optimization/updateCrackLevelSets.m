function [updated,Damage,info] = updateCrackLevelSets(LS,R,mask,parameters,baseDamage)
%UPDATECRACKLEVELSETS One energy-decreasing shape step at frozen displacement.
% The derivative of the discretized residual+crack energy is pulled back
% through the nearest-polyline distance and the exponential damage profile.
% Interior nodes move normally; both tip nodes may move in 2D. Backtracking
% accepts only a lower total energy and a noncrossing curve inside the ROI.
% All THREE level sets are then reconstructed as distance fields. This is
% a geometry-constrained shape optimization / reinitialization approach,
% not the independent damage PDE nor Rethore's Godunov propagation scheme.
p=phaseFieldParameters(parameters);
if nargin<5, baseDamage=zeros(size(R)); end
if ~isequal(size(R),LS.imageSize) || ~crackPathInDomain(LS.pathXY,mask)
    error('updateCrackLevelSets:Path','Initial path must lie inside the residual ROI.');
end
Damage=createDamageField(LS,mask,p,baseDamage);
[E,gd]=phaseFieldEnergy2D(Damage,R,mask,p);
[y,x]=find(mask);
[distance,segment,t,offset]=crackProjection2D(LS.pathXY,[x y]);
den=sqrt(distance.^2+p.distanceSmoothing^2);
factor=gd(mask).*sqrt(p.mu)/p.lc.*Damage(mask)./den;
n=size(LS.pathXY,1);
gradient=zeros(n,2);
for axis=1:2
    contribution=factor.*offset(:,axis);
    gradient(:,axis)=accumarray(segment,(1-t).*contribution,[n 1])+ ...
        accumarray(segment+1,t.*contribution,[n 1]);
end
% Remove interior tangential redistribution; endpoints retain tip motion.
for k=2:n-1
    tangent=LS.pathXY(k+1,:)-LS.pathXY(k-1,:);
    tangent=tangent/max(norm(tangent),eps);
    gradient(k,:)=gradient(k,:)-dot(gradient(k,:),tangent)*tangent;
end
maxGradient=max(vecnorm(gradient,2,2));
updated=LS;
info=struct('accepted',false,'energyBefore',E,'energyAfter',E, ...
    'step',0,'gradientNorm',maxGradient,'shapeGradient',gradient, ...
    'reason','stationary');
if maxGradient<1e-12, return; end
direction=-gradient/maxGradient;
info.reason='line-search-rejected';
for backtrack=0:p.maxBacktracks
    step=p.maxGeometryStep/2^backtrack;
    trialPath=LS.pathXY+step*direction;
    if ~crackPathInDomain(trialPath,mask), continue; end
    distances=crackProjection2D(trialPath,[x y]);
    trial=double(baseDamage);
    r=sqrt(distances.^2+p.distanceSmoothing^2)-p.distanceSmoothing;
    trial(mask)=exp(-sqrt(p.mu)*r/p.lc);
    trialEnergy=phaseFieldEnergy2D(trial,R,mask,p);
    if trialEnergy < E-1e-12*max(1,abs(E))
        updated=crackLevelSets2D(trialPath,size(R));
        Damage=trial;
        info.accepted=true; info.energyAfter=trialEnergy; info.step=step;
        info.reason='accepted'; return
    end
end
end
