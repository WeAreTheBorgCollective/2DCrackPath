function [BW, threshold, info] = valleyOtsuThreshold(rho, numBins)
%VALLEYOTSUTHRESHOLD Ng (2006), Eq. (8), for high-valued residual foreground.
% Maximizes (1-p(t))*(w1*mu1^2+w2*mu2^2), NOT weighted between-class variance.
% Nonnegative finite scores are scaled from [0,max(rho)] to integer levels
% 0:numBins-1 by rounding. Keeping zero fixed avoids an additive shift in
% Ng's objective. numBins (default 256) controls histogram quantization.
% threshold is the upper boundary of the selected bin in original GL units;
% BW uses rho >= threshold, exactly equivalent to quantized level > t.
% NaN/Inf are excluded. Empty/constant finite data return an empty mask.
% Reference: doi:10.1016/j.patrec.2006.03.009, Section 2.2.
if nargin < 2, numBins = 256; end
validateattributes(rho,{'numeric'},{'real','2d','nonempty'});
validateattributes(numBins,{'numeric'},{'scalar','integer','>=',2,'finite'});
valid = isfinite(rho);
v = double(rho(valid));
v = v(:);
if any(v < 0)
    error('valleyOtsuThreshold:Negative','Screening scores must be nonnegative.');
end
BW = false(size(rho));
threshold = NaN;
info = struct('probability',zeros(numBins,1),'objective',nan(numBins,1), ...
    'levels',(0:numBins-1).','thresholdLevel',NaN,'status','no-finite-data');
if isempty(v), return; end
threshold = max(v);
if min(v)==max(v)
    info.status = 'constant';
    return
end
scale = max(v)/(numBins-1);
q = min(numBins-1,max(0,round(v/scale)));
counts = accumarray(q+1,1,[numBins 1]);
p = counts/numel(v);
cumulativeCount = cumsum(counts);
w1 = cumulativeCount/numel(v); w2 = (numel(v)-cumulativeCount)/numel(v);
m1 = cumsum((0:numBins-1)'.*p); total = m1(end);
admissible = w1>0 & w2>0;
score = -inf(numBins,1);
score(admissible) = (1-p(admissible)).*( ...
    m1(admissible).^2./w1(admissible) + ...
    (total-m1(admissible)).^2./w2(admissible));
% First maximum is deterministic when multiple empty valley bins tie.
[~,index] = max(score);
threshold = (index-0.5)*scale;
BW(valid) = q > index-1;
info = struct('probability',p,'objective',score, ...
    'levels',(0:numBins-1)'*scale,'thresholdLevel',index-1,'status','ok');
end
