function tests = testScreenResidual2D
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'lib','optimization'));
testCase.TestData.root = root;
end

function testRectangularAgainstExplicitWindows(testCase)
R = reshape(sin(1:99),9,11);
R(5,6) = NaN; R(2,9) = Inf;
[rho,~,eligible] = screenResidual2D(R,[3 5],0,5,3);
expected = nan(size(R));
for y = 2:8
    for x = 3:9
        patch = R(y-1:y+1,x-2:x+2);
        if all(isfinite(patch(:)))
            expected(y,x) = mean(abs(patch(:)));
        end
    end
end
verifyEqual(testCase,isnan(rho),isnan(expected));
ok = isfinite(expected);
verifyEqual(testCase,rho(ok),expected(ok),'AbsTol',1e-12);
verifyEqual(testCase,eligible,ok & expected>0);
end

function testAbsoluteNotSquaredAndStrictThreshold(testCase)
[rho,c] = screenResidual2D([-3 3; 5 -5],[1 1],3,10,0);
verifyEqual(testCase,rho,[3 3;5 5]);
verifyEqual(testCase,c.meanAbsoluteResidual,[5;5]);
end

function testRankingAndSeparation(testCase)
R = zeros(12);
R(3,3)=10; R(3,4)=9; R(3,8)=8; R(8,3)=7;
[~,c] = screenResidual2D(R,[1 1],0,3,5);
verifyEqual(testCase,[c.x c.y],[3 3;8 3;3 8]);
verifyEqual(testCase,c.meanAbsoluteResidual,[10;8;7]);
end

function testEmptyAndInvalidInputs(testCase)
[~,c] = screenResidual2D(zeros(7),[3 5],0,5,5);
verifyEqual(testCase,height(c),0);
[rho,c] = screenResidual2D(nan(7),[3 5],0,5,5);
verifyTrue(testCase,all(isnan(rho(:))));
verifyEqual(testCase,height(c),0);
[~,c] = screenResidual2D(ones(7),[3 5],0,0,5);
verifyEqual(testCase,height(c),0);
verifyError(testCase,@()screenResidual2D(ones(7),[2 3],0,5,5), ...
    'screenResidual2D:WindowSize');
verifyError(testCase,@()screenResidual2D(ones(7),[9 3],0,5,5), ...
    'screenResidual2D:WindowSize');
end
