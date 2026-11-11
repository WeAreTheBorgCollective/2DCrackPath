function R = evaluateUndamagedResidual(kernel,f,g,a,measurementMask)
%EVALUATEUNDAMAGEDRESIDUAL Evaluate f(X)-g(X+u) on the ORIGINAL images.
% This interface intentionally has no damage argument.
R = kernel(single(f),single(g),reshape(a.',[],1),'R');
if ~isequal(size(R),size(f),size(g),size(measurementMask)) || ...
        any(~isfinite(a(:))) || any(~isfinite(R(measurementMask)))
    error('evaluateUndamagedResidual:Coverage', ...
        'Displacement must provide finite original-image residuals over the fixed ROI.');
end
end
