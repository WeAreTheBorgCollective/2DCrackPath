%% Zhu and Lu: independent damage PDE and weighted DIC, Section 2.7
% This is separate from the older prescribed-crack/level-set experiments.
% Initialize damage from extractResidualCrackPath; displacement starts at zero.
% Coordinates match the paper: g(x)-f(x-u(x)), on the DEFORMED image grid.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
p = paperPFParameters;
% Images are not Gaussian-filtered; use the original-image residual.
% Optional explicit calibration overrides, saved with every result.
if exist('paperParameterOverrides','var')
    overrideNames = fieldnames(paperParameterOverrides);
    for parameterIndex = 1:numel(overrideNames)
        name = overrideNames{parameterIndex};
        if ~isfield(p,name)
            error('paperPF:Parameter','Unknown parameter: %s',name);
        end
        p.(name) = paperParameterOverrides.(name);
    end
end
% This entry point always uses raw images, including if stale overrides exist.
p.displacementSmoothing = 0;
% Calibrate the Table 3 values for your experiment before quantitative use.
[fOriginal,gOriginal] = loadDICImages(projectRoot);
fOriginal = double(fOriginal);
gOriginal = double(gOriginal);
if ~isequal(size(fOriginal),size(gOriginal))
    error('paperPF:Images','Reference and deformed images must have the same size.');
end
grayMin = min([fOriginal(:);gOriginal(:)]);
grayRange = max([fOriginal(:);gOriginal(:)])-grayMin;
if grayRange<=0
    error('paperPF:Images','Images must have nonzero gray-level range.');
end
f = (fOriginal-grayMin)/grayRange;
g = (gOriginal-grayMin)/grayRange;
if p.invertImages
    f = 1-f;
    g = 1-g;
end

%% Rectangular specimen ROI, not a crack mask
% Reuse only the previous measurement rectangle to locate the specimen.
% Do not reuse the previous mesh exclusion or displacement coefficients.
initialGeometry = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'), ...
    'measurementMask');
[roiY,roiX] = find(initialGeometry.measurementMask);
measurementMask = false(size(f));
measurementMask(min(roiY):max(roiY),min(roiX):max(roiX)) = true;
model = paperPFModel(measurementMask,p.displacementSpacing);
plotBounds = [min(roiX) max(roiX) min(roiY) max(roiY)];
%% Extracted crack initializes damage on the deformed grid
crackFile = fullfile(projectRoot,'CrackPaths','CrackPathInit.mat');
crackData = load(crackFile);
[d,damageInitialization] = paperPFInitialDamage(model,crackData,p);
points = damageInitialization.pathXY;
plotBounds = [max(1,min(points(:,1))-20) min(size(f,2),max(points(:,1))+20) ...
    max(1,min(points(:,2))-20) min(size(f,1),max(points(:,2))+20)];
% Initial damage is free to change in subsequent independent PDE solves.
U = zeros(size(model.nodes));
r = g(model.ids)-f(model.ids);
Phi = zeros(size(r));
history = struct([]);
converged = false;
status = 'maximum-iterations';
outputRoot = fullfile(projectRoot,'PaperPhaseField');
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end
runFolder = fullfile(outputRoot,['Run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
mkdir(runFolder);
save(fullfile(runFolder,'Model.mat'),'model','p','grayMin','grayRange','plotBounds','damageInitialization','-v7.3');
initialDamage = zeros(model.imageSize);
initialDamage(model.ids) = d;
save(fullfile(runFolder,'InitialDamage.mat'),'initialDamage', ...
    'damageInitialization','measurementMask','p','-v7');
savePaperPFIteration(runFolder,0,model,U,d,r,Phi,history,p,gOriginal,plotBounds);
fprintf('Extracted-crack damage initialization; original images, no Gaussian filtering.\n');
firstDisplacement = [];
for iteration = 1:p.maxIterations
    % Displacement increment with FIXED damage and genuine quadrature weights.
    [U,deltaU,r,displacementInfo] = paperPFDisplacementSolve( ...
        model,f,g,U,d,p,iteration);
    if iteration==1
        firstDisplacement = model.N*U;
    end

    % Updated unweighted residual and strain drive an INDEPENDENT damage PDE.
    [Phi,drivingDetails] = paperPFDrivingForce(model,U,r,f,p);
    previousDamage = d;
    [d,damageInfo] = paperPFDamageSolve(model.K,model.mass,r.^2+Phi,d,p,model.colors);
    M = (1-d).^2;
    incrementPixels = model.N*deltaU;
    numerator = sqrt(sum(model.mass.*sum((M.*incrementPixels).^2,2)));
    denominator = sqrt(sum(model.mass.*sum((M.*firstDisplacement).^2,2)));
    if denominator>eps
        ratio = numerator/denominator;
    elseif numerator<=eps
        ratio = 0;
    else
        ratio = Inf;
    end
    entry = struct('iteration',iteration,'weightedIncrementRatio',ratio, ...
        'weightedIncrementNorm',numerator,'weightedFirstNorm',denominator, ...
        'maxDamage',max(d),'damageChange',norm(d-previousDamage,Inf), ...
        'damagePixels',nnz(d>=p.damageThreshold), ...
        'rawRMSE',sqrt(mean(r.^2)),'indicatorMax',max(drivingDetails.indicator), ...
        'displacement',displacementInfo,'damage',damageInfo);
    if iteration==1
        history = entry;
    else
        history(iteration) = entry;
    end
    savePaperPFIteration(runFolder,iteration,model,U,d,r,Phi,history,p,gOriginal,plotBounds);
    fprintf('Paper PF %02d: ratio %.4g, max d %.4f, damaged pixels %d, damage KKT %.3g\n', ...
        iteration,ratio,max(d),nnz(d>=p.damageThreshold),damageInfo.projectedKKT);
    if ~damageInfo.converged
        status = 'damage-solver-not-converged';
        break
    end
    if ~displacementInfo.accepted && displacementInfo.gradientNorm>1e-8
        status = 'displacement-line-search-stalled';
        break
    end
    if ratio<=p.tolerance
        converged = true;
        status = 'converged';
        break
    end
end

%% Save named final outputs and preserve all iterations in the timestamped run
finalFile = fullfile(runFolder,sprintf('Iteration_%02d.mat',iteration));
copyfile(finalFile,fullfile(runFolder,'Damage.mat'));
copyfile(finalFile,fullfile(outputRoot,'Damage.mat'));
save(fullfile(runFolder,'Damage.mat'),'converged','status','runFolder','-append');
save(fullfile(outputRoot,'Damage.mat'),'converged','status','runFolder','-append');
save(fullfile(runFolder,'Convergence.mat'),'history','converged','status','p','-v7');
save(fullfile(outputRoot,'LatestRun.mat'),'runFolder','converged','status','p','-v7');
fprintf('Paper PF-DIC status: %s. Outputs: %s\n',status,runFolder);
if ~converged
    warning('paperPF:NotConverged','Solver stopped with status: %s',status);
end
