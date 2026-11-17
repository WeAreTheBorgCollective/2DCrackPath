function [U,crackPathXY] = evaluateXDICDisplacement(result)
%EVALUATEXDICDISPLACEMENT Reconstruct the final enriched field without a solve.
prep = result.stat.prep;
crackPathXY = result.levelSets.pathXY;
assert(isequal(prep.crackPathXY,crackPathXY), ...
    'Final crack geometry and saved displacement basis do not match.');
assert(isequal(prep.enrichedNodes(:),result.solution.enrichedNodes(:)), ...
    'Enrichment ordering does not match the saved coefficients.');
coefficients = result.solution.coefficients;
assert(size(coefficients,1)==size(prep.basis,2) && size(coefficients,2)==2, ...
    'Incompatible displacement coefficient dimensions.');
pixels = prep.basis*coefficients;
imageSize = double(prep.imageSize(:)');
U = nan([imageSize 2]);
for component = 1:2
    field = nan(imageSize);
    field(prep.pixelIndex) = pixels(:,component);
    U(:,:,component) = field;
end
end
