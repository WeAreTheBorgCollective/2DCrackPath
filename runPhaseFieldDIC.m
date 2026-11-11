%% Phase-field DIC followed by screened-residual-gradient level-set evolution
% Run example2Dscreening, screenResidualImage, extractResidualCrackPath first.
% Configure Correli before running this script.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
if exist('correli_pathdef','file') == 2
    addpath(correli_pathdef);
end

%% Independent level-set settings
levelSetOptions = struct( ...
    'screeningWindowSize',[23 23], ...
    'velocityMethod','gradient', ...
    'gradientGain',1, ...       % px^2 / GL per unit virtual time
    'maxNormalStep',1, ...      % speed cap [px / virtual time]
    'duration',1, ...
    'smoothLength',3);

initial = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
crack = load(fullfile(projectRoot,'CrackPaths','CrackPathInit.mat'),'crackPathXY');
f = single(imread(fullfile(projectRoot,'Slike','reference.tif')));
g = single(imread(fullfile(projectRoot,'Slike','deformed.tif')));
if ~isequal(size(f),size(g),initial.imageSize)
    error('runPhaseFieldDIC:Size','Initialization and images must have the same size.');
end

measurementMask = initial.measurementMask;
Mesh = initial.Mesh;
a = initial.a;
parameters = phaseFieldParameters(struct( ...
    'grayScale',initial.grayScale,'lc',5,'mu',0.1,'outerIterations',3));
if ~crackPathInDomain(crack.crackPathXY,measurementMask)
    error('runPhaseFieldDIC:Path','Crack must be a noncrossing polyline inside the ROI.');
end

%% Construct all three level sets from the extracted crack
levelSets = crackLevelSets2D(crack.crackPathXY,size(f));
initialLevelSets = levelSets;
baseDamage = zeros(size(f));
Damage = createDamageField(levelSets,measurementMask,parameters,baseDamage);
kernel = kernel_corr('cubic',kernel_prep(Mesh,size(f)));

% Always evaluate residuals on original images. No damage weighting here.
R = evaluateUndamagedResidual(kernel,f,g,a,measurementMask);
damageFolder = fullfile(projectRoot,'Damage');
if ~isfolder(damageFolder)
    mkdir(damageFolder);
end
runStamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
evolutionFolder = fullfile(damageFolder,['Evolution_' runStamp]);
mkdir(evolutionFolder);
% Preserve the named outputs before this run replaces them.
for previousFile = {'Damage.mat','PhaseFieldResult.mat'}
    sourceFile = fullfile(damageFolder,previousFile{1});
    if isfile(sourceFile)
        copyfile(sourceFile,fullfile(evolutionFolder,['Previous_' previousFile{1}]));
    end
end
saveLevelSetStep(evolutionFolder,0,levelSets,Damage,a,R, ...
    measurementMask,levelSetOptions,'Initial conventional-DIC displacement');

energyHistory = phaseFieldEnergy2D(Damage,R,measurementMask,parameters);
updateHistory = cell(parameters.outerIterations,2);
for iteration = 1:parameters.outerIterations
    % 1. Solve PF-DIC for the CURRENT damage and retain the accepted a.
    [a,~,dicInfo] = correlatePhaseField( ...
        f,g,Mesh,a,R,Damage,measurementMask,parameters);

    % 2. Warp original g with that a and compare with original f.
    R = evaluateUndamagedResidual(kernel,f,g,a,measurementMask);

    % 3. Screen R, ascend its normal gradient, then rebuild phi and both tips.
    [levelSets,Damage,shapeInfo] = updateCrackLevelSetsRethore( ...
        levelSets,R,measurementMask,parameters,levelSetOptions,baseDamage);
    updateHistory(iteration,:) = {shapeInfo,dicInfo};
    energyHistory(end+1,1) = phaseFieldEnergy2D( ...
        Damage,R,measurementMask,parameters); %#ok<SAGROW>

    saveLevelSetStep(evolutionFolder,iteration,levelSets,Damage,a,R, ...
        measurementMask,levelSetOptions,'After shape update; residual from preceding PF solve');
    fprintf('Step %d: DIC accepted=%d; shape accepted=%d; screen %.5f -> %.5f GL\n', ...
        iteration,dicInfo.accepted,shapeInfo.accepted, ...
        shapeInfo.meanScreeningBefore,shapeInfo.meanScreeningAfter);
end

%% Synchronize the final displacement with the last damage geometry
[a,~,finalDICInfo] = correlatePhaseField( ...
    f,g,Mesh,a,R,Damage,measurementMask,parameters);
R = evaluateUndamagedResidual(kernel,f,g,a,measurementMask);
energyHistory(end+1,1) = phaseFieldEnergy2D(Damage,R,measurementMask,parameters);
saveLevelSetStep(evolutionFolder,parameters.outerIterations+1, ...
    levelSets,Damage,a,R,measurementMask,levelSetOptions,'Final PF solve');

save(fullfile(damageFolder,'Damage.mat'), ...
    'Damage','measurementMask','levelSets','initialLevelSets','parameters', ...
    'energyHistory','updateHistory','finalDICInfo','levelSetOptions','evolutionFolder','-v7');
save(fullfile(damageFolder,'PhaseFieldResult.mat'), ...
    'a','R','Mesh','parameters','energyHistory','finalDICInfo','evolutionFolder','-v7');

% Uses original images and refreshes InitialValues/DisplacementInitPF.mat.
comparePhaseFieldDIC;
fprintf('Level-set snapshots and figures: %s\n',evolutionFolder);
