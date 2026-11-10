function Damage = createDamageField(LS,measurementMask,parameters,baseDamage,outputFile)
%CREATEDAMAGEFIELD Geometry-constrained diffuse crack; optionally save a MAT.
% d=exp(-sqrt(mu)*(sqrt(distance^2+eps^2)-eps)/lc) inside the ROI.
% With eps -> 0 this is the exact straight-crack homogeneous profile for
% (mu/lc)d-lc*Delta(d)=0, with d=1 on the crack. Curved finite cracks use
% this as a distance-profile ansatz, NOT an exact solution of the paper's
% coupled damage PDE. eps regularizes the shape derivative at distance zero.
% Outside measurementMask, baseDamage is preserved (default zero). Image
% intensity data belong in f/g, never in the damage array. Those image pixels
% are left untouched by applyDamageToImages.
if nargin<3, parameters=struct(); end
p=phaseFieldParameters(parameters);
if ~islogical(measurementMask) || ~isequal(size(measurementMask),LS.imageSize)
    error('createDamageField:Mask','Mask must be logical and match the image size.');
end
if nargin<4 || isempty(baseDamage), baseDamage=zeros(LS.imageSize); end
validateattributes(baseDamage,{'numeric'},{'real','size',LS.imageSize});
Damage=double(baseDamage);
r=sqrt(LS.distance(measurementMask).^2+p.distanceSmoothing^2)-p.distanceSmoothing;
Damage(measurementMask)=exp(-sqrt(p.mu)*r/p.lc);
if nargin>=5 && ~isempty(outputFile)
    folder=fileparts(outputFile);
    if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
    levelSets=LS; parameters=p; %#ok<NASGU>
    save(outputFile,'Damage','measurementMask','levelSets','parameters','-v7');
end
end
