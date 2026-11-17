function score = scoreCrackPathInROI(path,image,mask)
% Mean screened score along uniformly spaced reference-curve samples.
s = [0;cumsum(vecnorm(diff(path),2,2))];
q = interp1(s,path,linspace(0,s(end),max(2,ceil(s(end)*2)+1))');
values = sampleROIImage(abs(double(image)),mask,q);
if any(~isfinite(values))
    score = -Inf;
else
    score = mean(values);
end
end
