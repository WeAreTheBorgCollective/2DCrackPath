function [Phi,details] = paperPFDrivingForce(model,U,r,f,p)
%PAPERPFDRIVINGFORCE Eqs.16,34, full strain tensor (no tensile-only split).
exx = model.Bx*U(:,1);
eyy = model.By*U(:,2);
exy = 0.5*(model.By*U(:,1)+model.Bx*U(:,2));
strainSquared = exx.^2+eyy.^2+2*exy.^2;
rMean = paperPFBlockMean(r,model.mask,p.residualBlock);
fMean = paperPFBlockMean(f(model.ids),model.mask,p.grayBlock);
indicator = fMean.*rMean.^2.*strainSquared;
Phi = p.w1*tanh((indicator-p.a1)/p.c1)+ ...
    p.w2*tanh((fMean-p.a2)/p.c2);
details = struct('strainSquared',strainSquared,'meanResidual',rMean, ...
    'meanGray',fMean,'indicator',indicator);
end
