function [energy,gradientD,parts] = phaseFieldEnergy2D(Damage,R,mask,parameters)
%PHASEFIELDENERGY2D Eq.8 of Zhu & Lu with Phi=0, unit pixel spacing.
% E=1/2 sum_ROI((1-d)^2*(R/grayScale)^2 + mu/lc*d^2)
%   +lc/2 sum_ROI-neighbor-bonds (d_j-d_i)^2.
% Uses only bonds whose BOTH ends lie inside the fixed ROI: natural no-flux
% boundaries. No artificial jumps to outside pixels contribute to energy.
p=phaseFieldParameters(parameters);
if ~islogical(mask) || ~isequal(size(Damage),size(R),size(mask)) || ~any(mask(:))
    error('phaseFieldEnergy2D:Size','Nonempty logical ROI and equal image sizes required.');
end
if any(~isfinite(R(mask))) || any(~isfinite(Damage(mask))) || ...
        any(Damage(mask)<0 | Damage(mask)>1)
    error('phaseFieldEnergy2D:Data','Finite residuals and damage in [0,1] are required inside ROI.');
end
d=zeros(size(Damage)); d(mask)=double(Damage(mask));
r=zeros(size(R)); r(mask)=double(R(mask))/p.grayScale;
dic=0.5*sum((1-d(mask)).^2.*r(mask).^2);
local=0.5*p.mu/p.lc*sum(d(mask).^2);
gradientD=zeros(size(d));
gradientD(mask)=-(1-d(mask)).*r(mask).^2+p.mu/p.lc*d(mask);
vertical=mask(1:end-1,:) & mask(2:end,:);
horizontal=mask(:,1:end-1) & mask(:,2:end);
dv=diff(d,1,1); dv(~vertical)=0;
dh=diff(d,1,2); dh(~horizontal)=0;
surface=0.5*p.lc*(sum(dv(:).^2)+sum(dh(:).^2));
gradientD(1:end-1,:)=gradientD(1:end-1,:)-p.lc*dv;
gradientD(2:end,:)=gradientD(2:end,:)+p.lc*dv;
gradientD(:,1:end-1)=gradientD(:,1:end-1)-p.lc*dh;
gradientD(:,2:end)=gradientD(:,2:end)+p.lc*dh;
energy=dic+local+surface;
parts=struct('dic',dic,'crack',local+surface,'localDamage',local,'gradientDamage',surface);
end
