function testSauvolaResidual
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setupPaths;
[x,y] = meshgrid(1:81);
scores = 10+.1*x;
scores(abs(y-41)<=2) = scores(abs(y-41)<=2)+100;
o = struct('windowSize',[21 21],'scalePercentile',100);
[b,t,info] = sauvolaResidualThreshold(scores,o);
assert(all(b(41,20:60)));
assert(~any(b(10,20:60)));
assert(isequal(b,sauvolaResidualThreshold(3*scores,o)));
% Independent window computation at one center, with explicit bright/dark polarity.
I = 1-scores/info.scaleGL;
window = I(31:51,31:51);
T = mean(window(:))*(1+.2*(std(window(:),1)/.5-1));
assert(abs(t(41,41)-info.scaleGL*(1-T))<1e-10);
assert(b(41,41)==(I(41,41)<T));
scores(35:37,40:42) = NaN;
b = sauvolaResidualThreshold(scores,o);
assert(~any(b(isnan(scores))));
assert(~any(sauvolaResidualThreshold(ones(40)),'all'));
assert(~any(sauvolaResidualThreshold(nan(40)),'all'));
fprintf('PASS: bright ridge, background, GL scale invariance, formula, NaNs and constants.\n');
end
