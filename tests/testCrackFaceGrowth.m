function testCrackFaceGrowth
path = [(10:20)' zeros(11,1)];
F = repmat(eye(2),1,1,2);
% Right-hand rays meet at (25,0); left-hand rays diverge.
minus = [path(:,1) -(25-path(:,1))/5];
plus = [path(:,1) (25-path(:,1))/5];
[q,info] = growCrackFromFaceTraces(path,minus,plus,F);
assert(~info(1).accepted && info(2).accepted);
assert(norm(q(end,:)-[25 0])<1e-10);
assert(abs(info(2).extension-5)<1e-10);
% Far intersection is capped to 10 reference pixels.
minus(:,2) = -(50-path(:,1))/5;
plus(:,2) = (50-path(:,1))/5;
[q,info] = growCrackFromFaceTraces(path,minus,plus,F);
assert(norm(q(end,:)-[30 0])<1e-10);
assert(abs(info(2).uncappedExtension-30)<1e-10);
% Translation and affine stretch of the deformed image must not scale growth.
A = [2 0.3;0 1.5];
offset = [70 -20];
minus = [path(:,1) -(25-path(:,1))/5]*A'+offset;
plus = [path(:,1) (25-path(:,1))/5]*A'+offset;
[q,~] = growCrackFromFaceTraces(path,minus,plus,repmat(A,1,1,2));
assert(norm(q(end,:)-[25 0])<1e-10);
% Reverse ordering preserves geometry and exchanges the growing tip.
[qReverse,~] = growCrackFromFaceTraces(flipud(path),flipud(minus), ...
    flipud(plus),repmat(A,1,1,2));
assert(max(abs(qReverse(:)-reshape(flipud(q),[],1)))<1e-10);
% Parallel traces and a coincident tip never grow.
minus = path+[0 -1];
plus = path+[0 1];
[q,info] = growCrackFromFaceTraces(path,minus,plus,F);
assert(isequal(q,path) && ~any([info.accepted]));
minus(:,2) = -(20-path(:,1))/5;
plus(:,2) = (20-path(:,1))/5;
[q,info] = growCrackFromFaceTraces(path,minus,plus,F);
assert(isequal(q,path) && ~any([info.accepted]));
fprintf('Crack-face intersection, cap and coordinate-frame tests passed.\n');
end
