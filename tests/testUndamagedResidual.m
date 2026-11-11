function tests = testUndamagedResidual
tests = functiontests(localfunctions);
end

function setupOnce(~)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'lib','corr'));
end

function testOriginalImagesAndDisplacementPacking(tc)
f = single(reshape(1:30,5,6));
g = f+10;
a = [1 2 3;4 5 6];
mask = true(size(f));
kernel = @(reference,deformed,packed,mode) checkInputs( ...
    reference,deformed,packed,mode,f,g,a);
R = evaluateUndamagedResidual(kernel,f,g,a,mask);
verifyEqual(tc,R,f-g);
end

function R = checkInputs(reference,deformed,packed,mode,f,g,a)
assert(isequal(reference,f));
assert(isequal(deformed,g));
assert(isequal(packed,reshape(a.',[],1)));
assert(strcmp(mode,'R'));
R = reference-deformed;
end

function testCoverageLossRejected(tc)
f = zeros(5,6,'single');
kernel = @(varargin) nan(size(f));
verifyError(tc,@()evaluateUndamagedResidual(kernel,f,f,zeros(2,3),true(size(f))), ...
    'evaluateUndamagedResidual:Coverage');
end
