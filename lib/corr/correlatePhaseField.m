function [a,R,info] = correlatePhaseField(f,g,Mesh,initialA,initialR,Damage,mask,parameters)
%CORRELATEPHASEFIELD Use existing correlate on damage-weighted images.
% Re-evaluate candidate displacements on ORIGINAL images using the kernel
% handle returned by correlate. Accept only if the actual mechanical-free
% PF energy decreases on the SAME fixed ROI. Reject any ROI coverage loss.
% BCC and mechanical regularization are disabled. Images use original GL;
% normalization occurs only in phaseFieldEnergy2D via grayScale.
p=phaseFieldParameters(parameters);
if exist('correlate','file')~=2
    error('correlatePhaseField:Dependency','Configure Correli (addpath(correli_pathdef)) first.');
end
[fPF,gPF]=applyDamageToImages(single(f),single(g),Damage,mask);
[candidate,~,~,stat]=correlate(fPF,gPF,Mesh,initialA, ...
    'plot_freq',0,'verbose_freq',0,'max_iter',100,'conv_limit',1e-6, ...
    'bcc',false,'mechreg',[]);
oldEnergy=phaseFieldEnergy2D(Damage,initialR,mask,p);
rawR=stat.corr_func(single(f),single(g),reshape(candidate.',[],1),'R');
newEnergy=Inf;
if all(isfinite(rawR(mask))) && all(isfinite(candidate(:)))
    newEnergy=phaseFieldEnergy2D(Damage,rawR,mask,p);
end
a=initialA; R=initialR;
info=struct('accepted',newEnergy<=oldEnergy,'energyBefore',oldEnergy, ...
    'candidateEnergy',newEnergy,'energyAfter',oldEnergy,'solverStatus',stat.conv_stat);
if info.accepted
    a=candidate; R=rawR; info.energyAfter=newEnergy;
end
end
