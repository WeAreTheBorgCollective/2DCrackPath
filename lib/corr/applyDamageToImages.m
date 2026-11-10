function [fPF,gPF] = applyDamageToImages(f,g,Damage,mask)
%APPLYDAMAGETOIMAGES Apply sqrt(M)=1-d to BOTH images, inside the ROI only.
% For f(X)-g(X+u), fixed-grid masking gives (1-d(X))*f(X) -
% (1-d(X+u))*g(X+u), not exactly (1-d(X))*(f(X)-g(X+u)).
% This requested image-only surrogate is exact at u=0 (or constant damage)
% but is NOT a replacement for a weighted correlator. Original-image
% residuals must be used for actual PF-energy evaluation/acceptance.
% Inputs are single/double so that outside-ROI pixel values are unchanged.
validateattributes(f,{'single','double'},{'real','2d'});
validateattributes(g,{'single','double'},{'real','size',size(f)});
if ~islogical(mask) || ~isequal(size(f),size(Damage),size(mask))
    error('applyDamageToImages:Size','Images, damage and logical ROI must match.');
end
if any(~isfinite(Damage(mask))) || any(Damage(mask)<0 | Damage(mask)>1)
    error('applyDamageToImages:Damage','ROI damage must be finite and in [0,1].');
end
fPF=f; gPF=g;
w=1-Damage(mask);
fPF(mask)=f(mask).*cast(w,'like',f);
gPF(mask)=g(mask).*cast(w,'like',g);
end
