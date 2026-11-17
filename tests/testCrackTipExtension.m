function testCrackTipExtension
% Endpoint extension obeys forward direction, 1 px tolerance and component.
binary = false(40,60);
binary(18:22,10:40) = true;
binary(18:22,48:55) = true; % disconnected distractor on the same ray
path = [(15:35)' repmat(20,21,1)];
[out,mask,info] = extendCrackPathToMask(path,binary,8);
assert(isequal(out(1,:),[10 20]) && isequal(out(end,:),[40 20]));
assert(isequal(info.extensionLengths,[5;5]));
assert(all(mask(20,10:40)));
assert(all(vecnorm(diff(out),2,2)>0));
% Beyond the endpoint, same-component pixels too far from the ray are ignored.
binary(22,40:45) = true;
[out,~,~] = extendCrackPathToMask(path,binary,8);
assert(isequal(out(end,:),[40 20]));
% A pixel at exactly 1 px perpendicular distance is eligible.
binary(21,41) = true;
[out,~,~] = extendCrackPathToMask(path,binary,8);
assert(isequal(out(end,:),[41 21]));
% Reversing path order reverses the same geometric extension.
[reverse,~,~] = extendCrackPathToMask(flipud(path),binary,8);
assert(max(abs(reverse(:)-reshape(flipud(out),[],1)))<1e-12);
% Diagonal ray distances are Euclidean, not row/column offsets.
binary = false(40);
for k = 5:30
    binary(k,k) = true;
end
[out,~,~] = extendCrackPathToMask([(10:25)' (10:25)'],binary,6);
assert(isequal(out(1,:),[5 5]) && isequal(out(end,:),[30 30]));
% A path already reaching its region ends gains no duplicate vertices.
path = [(5:30)' (5:30)'];
[out,~,info] = extendCrackPathToMask(path,binary,6);
assert(isequal(out,path) && ~any(info.extensionLengths));
fprintf('Crack tip extension tests passed.\n');
end
