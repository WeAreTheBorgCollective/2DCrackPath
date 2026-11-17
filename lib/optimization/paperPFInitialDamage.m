function [d,info] = paperPFInitialDamage(model,crackData,p)
%PAPERPFINITIALDAMAGE Seed independent damage from the extracted finite crack.
% PF uses deformed image coordinates. Never interpret reference XY as deformed.
% Homogeneous straight-crack profile exp(-sqrt(mu)*distance/lc), used only
% as an initial guess, not a prescribed field or an irreversibility constraint.
if ~isfield(crackData,'crackPathDeformedXY')
    error('paperPF:CrackFrame', ...
        'CrackPathInit.mat needs crackPathDeformedXY. Run extractResidualCrackPath first.');
end
path = double(crackData.crackPathDeformedXY);
validateattributes(path,{'double'},{'real','finite','2d','ncols',2,'nonempty'});
if size(path,1)<2 || any(vecnorm(diff(path),2,2)<=eps)
    error('paperPF:CrackPath','Expected an ordered nondegenerate crack polyline.');
end
if isfield(crackData,'imageSize') && ~isequal(double(crackData.imageSize(:)'),double(model.imageSize(:)'))
    error('paperPF:CrackSize','Extracted crack and PF image sizes differ.');
end
validateattributes(p.mu,{'numeric'},{'scalar','finite','positive'});
validateattributes(p.lc,{'numeric'},{'scalar','finite','positive'});
if ~crackPathInDomain(path,model.mask)
    error('paperPF:CrackDomain','Extracted deformed crack must lie in the PF measurement domain.');
end
distance = crackProjection2D(path,[model.x model.y]);
d = exp(-sqrt(p.mu)*distance/p.lc);
info = struct('source','CrackPaths/CrackPathInit.mat', ...
    'pathField','crackPathDeformedXY','coordinateFrame','deformed', ...
    'pathXY',path,'profile','exp(-sqrt(mu)*distance/lc)', ...
    'decayLengthPixels',p.lc/sqrt(p.mu), ...
    'initialOnly',true,'imageFiltering','none');
end
