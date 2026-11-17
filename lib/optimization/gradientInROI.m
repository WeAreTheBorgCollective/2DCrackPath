function [gx,gy] = gradientInROI(image,mask)
% Central differences where possible, one-sided differences at ROI borders.
% No difference ever uses a pixel outside the finite ROI.
a = double(image);
valid = mask & isfinite(a);
a(~valid) = NaN;
left = [nan(size(a,1),1) a(:,1:end-1)];
right = [a(:,2:end) nan(size(a,1),1)];
above = [nan(1,size(a,2));a(1:end-1,:)];
below = [a(2:end,:);nan(1,size(a,2))];
gx = derivative(a,left,right,valid);
gy = derivative(a,above,below,valid);
end

function d = derivative(a,minus,plus,valid)
d = nan(size(a));
d(valid) = 0;
both = valid & isfinite(minus) & isfinite(plus);
forward = valid & ~isfinite(minus) & isfinite(plus);
backward = valid & isfinite(minus) & ~isfinite(plus);
d(both) = (plus(both)-minus(both))/2;
d(forward) = plus(forward)-a(forward);
d(backward) = a(backward)-minus(backward);
end
