%% Zhu and Lu: independent damage PDE and weighted DIC, Section 2.7
% This is separate from the older prescribed-crack/level-set experiments.
% No extracted crack or displacement is used to initialize the algorithm.
% Coordinates match the paper: g(x)-f(x-u(x)), on the DEFORMED image grid.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
p = paperPFParameters;
% First displacement solve: imgaussfilt(f,2), imgaussfilt(g,2), then originals.
% Filtering is registration continuation on the same grid, not downsampling.
% Damage and final residual use the unfiltered images after this refinement.
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
% Calibrate the Table 3 values for your experiment before quantitative use.
fOriginal = double(imread(fullfile(projectRoot,'Slike','reference.tif')));
gOriginal = double(imread(fullfile(projectRoot,'Slike','deformed.tif')));
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
% Do not reuse its crack exclusion, displacement, damage or extracted path.
initialGeometry = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'), ...
    'measurementMask');
[roiY,roiX] = find(initialGeometry.measurementMask);
measurementMask = false(size(f));
measurementMask(min(roiY):max(roiY),min(roiX):max(roiX)) = true;
model = paperPFModel(measurementMask,p.displacementSpacing);
plotBounds = [min(roiX) max(roiX) min(roiY) max(roiY)];
% Optional plot-only zoom. It never enters either numerical subproblem.
crackFile = fullfile(projectRoot,'CrackPaths','CrackPathInit.mat');
if isfile(crackFile)
    displayCrack = load(crackFile,'crackPathDeformedXY');
    if isfield(displayCrack,'crackPathDeformedXY') && ~isempty(displayCrack.crackPathDeformedXY)
        points = displayCrack.crackPathDeformedXY;
        plotBounds = [max(1,min(points(:,1))-20) min(size(f,2),max(points(:,1))+20) ...
            max(1,min(points(:,2))-20) min(size(f,1),max(points(:,2))+20)];
    end
end

%% Section 2.7 initialization: u0=0 and d0=0
U = zeros(size(model.nodes));
d = zeros(numel(model.ids),1);
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
save(fullfile(runFolder,'Model.mat'),'model','p','grayMin','grayRange','plotBounds','-v7.3');
savePaperPFIteration(runFolder,0,model,U,d,r,Phi,history,p,gOriginal,plotBounds);
fprintf('Gaussian displacement initialization scales [px]: %s\n',mat2str(p.displacementSmoothing));
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
