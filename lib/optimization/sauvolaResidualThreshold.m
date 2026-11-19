function [BW,thresholdGL,info] = sauvolaResidualThreshold(scores,options)
%SAUVOLARESIDUALTHRESHOLD Sauvola on inverted nonnegative residual scores.
% I=1-min(scores/scale,1); T=m*(1+k*(s/R-1)); foreground I<T.
% Bright residual equivalent: min(scores,scale)>scale*(1-T).
% scale is a recorded robust percentile; no smoothing of source images.
% Statistics use finite ROI pixels only, population standard deviation.
% Boundary windows need minCoverage of their full requested area.
% Sauvola & Pietikainen (2000), doi:10.1016/S0031-3203(99)00055-2.
if nargin<2
    options = struct;
end
o = struct('windowSize',[51 51],'k',0.2,'dynamicRange',0.5, ...
    'scalePercentile',99.5,'minCoverage',0.5);
names = fieldnames(options);
for j = 1:numel(names)
    if ~isfield(o,names{j})
        error('sauvola:Option','Unknown option %s',names{j});
    end
    o.(names{j}) = options.(names{j});
end
validateattributes(scores,{'numeric'},{'2d','real','nonempty'});
validateattributes(o.windowSize,{'numeric'},{'vector','numel',2,'integer','positive','finite'});
if any(mod(o.windowSize,2)~=1)
    error('sauvola:Window','Window dimensions must be odd.');
end
validateattributes(o.k,{'numeric'},{'scalar','finite','>',0,'<',1});
validateattributes(o.dynamicRange,{'numeric'},{'scalar','finite','positive'});
validateattributes(o.scalePercentile,{'numeric'},{'scalar','finite','>',0,'<=',100});
validateattributes(o.minCoverage,{'numeric'},{'scalar','finite','>',0,'<=',1});
valid = isfinite(scores);
values = double(scores(valid));
if any(values<0)
    error('sauvola:Negative','Use absolute residuals or nonnegative screening scores.');
end
BW = false(size(scores));
thresholdGL = nan(size(scores));
info = struct('options',o,'scaleGL',NaN,'status','no-finite-data', ...
    'validCenters',0,'clippedPixels',0);
if isempty(values)
    return
end
if min(values)==max(values)
    info.status = 'constant';
    return
end
scale = prctile(values,o.scalePercentile);
if scale<=0
    scale = max(values);
end
I = zeros(size(scores));
I(valid) = 1-min(double(scores(valid))/scale,1);
h = ones(o.windowSize(1),1);
w = ones(1,o.windowSize(2));
count = conv2(h,w,double(valid),'same');
average = conv2(h,w,I,'same')./max(count,1);
variance = max(0,conv2(h,w,I.^2,'same')./max(count,1)-average.^2);
T = average.*(1+o.k*(sqrt(variance)/o.dynamicRange-1));
eligible = valid & count>=o.minCoverage*prod(o.windowSize);
thresholdGL(eligible) = scale*(1-T(eligible));
BW(eligible) = I(eligible)<T(eligible);
info.scaleGL = scale;
info.status = 'ok';
info.validCenters = nnz(eligible);
info.clippedPixels = nnz(values>scale);
end
