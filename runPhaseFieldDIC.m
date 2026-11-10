%% Geometry-constrained phase-field DIC, without mechanical driving energy
% Workflow: example2Dscreening -> extractResidualCrackPath -> this script.
% Requires Correli on path. Generated MAT files are saved relative to project.
projectRoot=fileparts(mfilename('fullpath')); addpath(projectRoot); setupPaths;
if exist('correli_pathdef','file')==2, addpath(correli_pathdef); end
initial=load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
crack=load(fullfile(projectRoot,'CrackPaths','CrackPathInit.mat'),'crackPathXY');
f=single(imread(fullfile(projectRoot,'Slike','reference.tif')));
g=single(imread(fullfile(projectRoot,'Slike','deformed.tif')));
if ~isequal(size(f),size(g),initial.imageSize)
    error('runPhaseFieldDIC:Size','Initialization and current images must have the same size.');
end
measurementMask=initial.measurementMask;
parameters=phaseFieldParameters(struct('grayScale',initial.grayScale, ...
    'lc',5,'mu',0.1,'outerIterations',3,'maxGeometryStep',0.25));
a=initial.a; R=initial.R; Mesh=initial.Mesh;
if ~crackPathInDomain(crack.crackPathXY,measurementMask)
    error('runPhaseFieldDIC:Path','Crack must be a noncrossing open polyline inside the DIC ROI.');
end
levelSets=crackLevelSets2D(crack.crackPathXY,size(f));
initialLevelSets=levelSets;
baseDamage=zeros(size(f)); % Outside ROI remains zero; input images unchanged there.
damageFolder=fullfile(projectRoot,'Damage');
Damage=createDamageField(levelSets,measurementMask,parameters,baseDamage, ...
    fullfile(damageFolder,'Damage.mat'));
energyHistory=phaseFieldEnergy2D(Damage,R,measurementMask,parameters);
updateHistory=cell(parameters.outerIterations,2);
for iteration=1:parameters.outerIterations
    [levelSets,Damage,shapeInfo]=updateCrackLevelSets( ...
        levelSets,R,measurementMask,parameters,baseDamage);
    [a,R,dicInfo]=correlatePhaseField(f,g,Mesh,a,R,Damage,measurementMask,parameters);
    updateHistory(iteration,:)={shapeInfo,dicInfo};
    energyHistory(end+1,1)=phaseFieldEnergy2D(Damage,R,measurementMask,parameters); %#ok<SAGROW>
    fprintf('PF step %d: E %.9g -> %.9g; shape accepted=%d, DIC accepted=%d\n', ...
        iteration,energyHistory(end-1),energyHistory(end),shapeInfo.accepted,dicInfo.accepted);
end
save(fullfile(damageFolder,'Damage.mat'),'Damage','measurementMask','levelSets', ...
    'initialLevelSets','parameters','energyHistory','updateHistory','-v7');
save(fullfile(damageFolder,'PhaseFieldResult.mat'),'a','R','Mesh','parameters', ...
    'energyHistory','-v7');
figure('Color','w','Position',[60 60 1400 650]);
tiledlayout(2,3,'TileSpacing','compact');
ax=nexttile; imagesc(ax,levelSets.phi); axis(ax,'image'); colorbar(ax); title(ax,'Crack level set');
ax=nexttile; imagesc(ax,levelSets.psiStart); axis(ax,'image'); colorbar(ax); title(ax,'Start-tip level set');
ax=nexttile; imagesc(ax,levelSets.psiEnd); axis(ax,'image'); colorbar(ax); title(ax,'End-tip level set');
ax=nexttile; imagesc(ax,Damage,[0 1]); axis(ax,'image'); colorbar(ax); title(ax,'Damage d');
ax=nexttile; im=imagesc(ax,R); im.AlphaData=measurementMask; axis(ax,'image'); hold(ax,'on');
plot(ax,initialLevelSets.pathXY(:,1),initialLevelSets.pathXY(:,2),'w--');
plot(ax,levelSets.pathXY(:,1),levelSets.pathXY(:,2),'r-'); title(ax,'Initial / updated crack');
ax=nexttile; plot(ax,0:numel(energyHistory)-1,energyHistory,'o-'); grid(ax,'on');
xlabel(ax,'Staggered step'); ylabel(ax,'PF energy'); title(ax,'Accepted total energy');
