function testXDICDisplacementEvaluation
path = [1 1;2 2];
prep = struct('basis',sparse([1 0 -0.5;0 1 0.5]), ...
    'pixelIndex',[1;4],'imageSize',[2 2],'enrichedNodes',1,'crackPathXY',path);
solution = struct('coefficients',[2 3;4 5;6 8],'enrichedNodes',1);
result = struct('stat',struct('prep',prep), ...
    'levelSets',struct('pathXY',path),'solution',solution);
[U,p] = evaluateXDICDisplacement(result);
assert(isequal(p,path));
assert(isequal([U(1,1,1) U(1,1,2)],[-1 -1]));
assert(isequal([U(2,2,1) U(2,2,2)],[7 9]));
assert(all(isnan(U(1,2,:))));
result.levelSets.pathXY(2,1) = 3;
rejected = false;
try
    evaluateXDICDisplacement(result);
catch
    rejected = true;
end
assert(rejected,'Mismatched crack geometry must be rejected.');
fprintf('Enriched displacement reconstruction tests passed.\n');
end
