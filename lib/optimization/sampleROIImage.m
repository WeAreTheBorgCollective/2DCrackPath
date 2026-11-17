function values = sampleROIImage(image,mask,points)
% Bilinear interpolation renormalized over finite ROI neighbors.
% Query centers must themselves lie in the ROI by nearest-pixel membership,
% matching crackPathInDomain; missing neighbors never act as zero samples.
valid = mask & isfinite(image);
filled = double(image);
filled(~valid) = 0;
weight = interp2(double(valid),points(:,1),points(:,2),'linear',0);
numerator = interp2(filled,points(:,1),points(:,2),'linear',0);
inside = interp2(double(mask),points(:,1),points(:,2),'nearest',0)>0.5;
values = nan(size(points,1),1);
use = inside & weight>1e-12;
values(use) = numerator(use)./weight(use);
end
