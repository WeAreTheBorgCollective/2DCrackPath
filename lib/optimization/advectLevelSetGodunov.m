function [phi,info] = advectLevelSetGodunov(phi,F,duration,cfl)
%ADVECTLEVELSETGODUNOV Rethore et al. Eqs.26-32, unit pixel spacing.
% phi_t + F*|grad(phi)|=0; one-sided Godunov differences and explicit Euler.
% F is held fixed during this virtual-time interval. Conservative 2D CFL.
validateattributes(phi,{'numeric'},{'real','finite','2d','nonempty'});
validateattributes(F,{'numeric'},{'real','finite','size',size(phi)});
validateattributes(duration,{'numeric'},{'scalar','real','finite','nonnegative'});
if nargin<4, cfl=0.45; end
validateattributes(cfl,{'numeric'},{'scalar','real','positive','<=',0.5});
speed=max(abs(F(:))); n=max(1,ceil(2*duration*speed/cfl));
dt=duration/n; phi=double(phi);
for k=1:n
    dxm=phi-phi(:,[1 1:end-1]); dxp=phi(:,[2:end end])-phi;
    dym=phi-phi([1 1:end-1],:); dyp=phi([2:end end],:)-phi;
    plus=sqrt(max(max(dxm,0).^2,min(dxp,0).^2)+ ...
              max(max(dym,0).^2,min(dyp,0).^2));
    minus=sqrt(max(min(dxm,0).^2,max(dxp,0).^2)+ ...
               max(min(dym,0).^2,max(dyp,0).^2));
    phi=phi-dt*(max(F,0).*plus+min(F,0).*minus);
end
info=struct('substeps',n,'dt',dt,'duration',duration,'maxSpeed',speed);
end
