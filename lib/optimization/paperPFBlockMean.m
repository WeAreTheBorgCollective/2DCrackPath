function averaged = paperPFBlockMean(values,mask,blockSize)
%PAPERPFBLOCKMEAN Eq.34: non-overlapping blocks, NOT sliding |residual| windows.
% Blocks are anchored at the ROI bounding-box origin. Boundary blocks average
% their ROI pixels only. Signed residuals are averaged before squaring.
[y,x] = find(mask);
labels = floor((y-min(y))/blockSize)+ ...
    (ceil((max(y)-min(y)+1)/blockSize))*floor((x-min(x))/blockSize)+1;
sums = accumarray(labels,double(values));
counts = accumarray(labels,1);
averaged = sums(labels)./counts(labels);
end
