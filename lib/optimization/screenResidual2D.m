function [rho, candidates, eligible] = screenResidual2D(R, windowSize, threshold, maxCandidates, minDistance)
%SCREENRESIDUAL2D Rank centered rectangular windows by mean absolute residual.
%   windowSize = [height width] in pixels; both dimensions must be odd.
%   rho(y,x) = mean(abs(R(window centered at row y, column x))).
%   This is the screening statistic in Report, Section 4.2 (units: GL).
%   Only complete, finite windows are eligible; other scores are NaN.
%   Candidates exceed threshold strictly and are ranked by decreasing rho.
%   Greedy selection enforces Euclidean center separation >= minDistance.
%   Ties are resolved by MATLAB linear index (column first, then row).
%   Coordinates x,y are column,row pixel indices, starting at 1.
%   Candidates are alternative trial centers, not accepted cracks. No
%   PF-DIC re-solve, energy acceptance or level-set evolution is performed.

validateattributes(R, {'numeric'}, {'real','2d','nonempty'}, mfilename, 'R');
validateattributes(windowSize, {'numeric'}, ...
    {'real','finite','vector','numel',2,'integer','positive'}, mfilename, 'windowSize');
windowSize = double(windowSize(:).');
if any(mod(windowSize,2) ~= 1) || any(windowSize > size(R))
    error('screenResidual2D:WindowSize', ...
        'Window dimensions must be odd and no larger than the image.');
end
validateattributes(threshold, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(maxCandidates, {'numeric'}, {'real','finite','scalar','integer','nonnegative'});
validateattributes(minDistance, {'numeric'}, {'real','finite','scalar','nonnegative'});

R = double(R);
valid = isfinite(R);
values = abs(R);
values(~valid) = 0;
% Separable box sums avoid storing every extracted window.
h = windowSize(1); w = windowSize(2);
sums = conv2(ones(h,1), ones(1,w), values, 'same');
counts = conv2(ones(h,1), ones(1,w), double(valid), 'same');
rho = sums / (h*w);
rho(counts ~= h*w) = NaN;
eligible = isfinite(rho) & rho > threshold;

ids = find(eligible);
[~, order] = sortrows([-rho(ids), double(ids)], [1 2]);
ids = ids(order);
selected = zeros(min(double(maxCandidates), numel(ids)), 1);
nSelected = 0;
for k = 1:numel(ids)
    if nSelected == numel(selected)
        break
    end
    [y,x] = ind2sub(size(R), ids(k));
    [previousY,previousX] = ind2sub(size(R), selected(1:nSelected));
    if all(hypot(x-previousX, y-previousY) >= minDistance)
        nSelected = nSelected + 1;
        selected(nSelected) = ids(k);
    end
end
selected = selected(1:nSelected);
[y,x] = ind2sub(size(R), selected);
candidates = table((1:nSelected).', x, y, rho(selected), ...
    'VariableNames', {'rank','x','y','meanAbsoluteResidual'});
end
