%% Resume the FINAL accepted PF displacement and refine the crack by screened |R|
% Does not modify extractResidualCrackPath. Original outputs are backed up.
projectRoot=fileparts(mfilename('fullpath')); addpath(projectRoot); setupPaths;
if exist('correli_pathdef','file')==2, addpath(correli_pathdef); end
initial=load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
result=load(fullfile(projectRoot,'Damage','PhaseFieldResult.mat'));
saved=load(fullfile(projectRoot,'Damage','Damage.mat'));
f=single(imread(fullfile(projectRoot,'Slike','reference.tif')));
g=single(imread(fullfile(projectRoot,'Slike','deformed.tif')));
if ~isequal(size(f),size(g),initial.imageSize,saved.levelSets.imageSize) || ...
        ~isequal(initial.Mesh.pos,result.Mesh.pos) || ~isequal(initial.Mesh.con,result.Mesh.con)
    error('refineCrackPathRethore:Data','Images, ROI, level sets and displacement mesh must match.');
end
a=result.a; Mesh=result.Mesh; parameters=phaseFieldParameters(result.parameters);
measurementMask=saved.measurementMask; levelSets=saved.levelSets;
beforeLevelSets=levelSets; beforeA=a;
kernel=kernel_corr('cubic',kernel_prep(Mesh,size(f)));
R=kernel(f,g,reshape(a.',[],1),'R'); beforeResidual=R;
if any(~isfinite(R(measurementMask)))
    error('refineCrackPathRethore:Coverage','Final PF displacement loses original ROI coverage.');
end
options=struct('screeningWindowSize',[23 23], ... % independent of initial window
    'searchRadius',8,'sampleStep',0.25,'smoothLength',3,'duration',1);
iterations=3;
baseDamage=saved.Damage; Damage=createDamageField(levelSets,measurementMask,parameters,baseDamage);
energyHistory=phaseFieldEnergy2D(Damage,R,measurementMask,parameters);
updateHistory=cell(iterations,2);
for k=1:iterations
    % Recompute from current PF a on ORIGINAL f/g, never weighted images.
    R=kernel(f,g,reshape(a.',[],1),'R');
    [levelSets,Damage,shapeInfo]=updateCrackLevelSetsRethore( ...
        levelSets,R,measurementMask,parameters,options,baseDamage);
    [a,R,dicInfo]=correlatePhaseField(f,g,Mesh,a,R,Damage,measurementMask,parameters);
    updateHistory(k,:)={shapeInfo,dicInfo};
    energyHistory(end+1,1)=phaseFieldEnergy2D(Damage,R,measurementMask,parameters); %#ok<SAGROW>
    fprintf('Rethore step %d: mean screening score on crack %.4f -> %.4f GL; shape=%d; DIC=%d; E=%.6f\n', ...
        k,shapeInfo.meanScreeningBefore,shapeInfo.meanScreeningAfter,shapeInfo.accepted, ...
        dicInfo.accepted,energyHistory(end));
end
folder=fullfile(projectRoot,'Damage');
% Keep the starting solution for direct comparison, without replacing inputs.
backupFolder=fullfile(folder,['BeforeRethore_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
mkdir(backupFolder);
copyfile(fullfile(folder,'Damage.mat'),fullfile(backupFolder,'Damage.mat'));
copyfile(fullfile(folder,'PhaseFieldResult.mat'),fullfile(backupFolder,'PhaseFieldResult.mat'));
initialLevelSets=saved.initialLevelSets;
save(fullfile(folder,'Damage.mat'),'Damage','measurementMask','levelSets', ...
    'initialLevelSets','parameters','energyHistory','updateHistory','options','-v7');
save(fullfile(folder,'PhaseFieldResult.mat'),'a','R','Mesh','parameters','energyHistory','-v7');
save(fullfile(folder,'RethoreComparison.mat'),'beforeLevelSets','levelSets', ...
    'beforeA','a','beforeResidual','R','options','updateHistory','energyHistory','-v7');
figure('Color','w','Position',[50 50 1400 500]); tiledlayout(1,2,'TileSpacing','compact');
old=beforeLevelSets.pathXY; new=levelSets.pathXY;
bounds=[min([old;new],[],1)-15;max([old;new],[],1)+15];
for panel=1:2
    ax=nexttile; im=imagesc(ax,abs(double(beforeResidual)));
    im.AlphaData=measurementMask; axis(ax,'image'); hold(ax,'on');
    plot(ax,old(:,1),old(:,2),'w--','LineWidth',1.5);
    plot(ax,new(:,1),new(:,2),'r-','LineWidth',1.5);
    colorbar(ax); xlabel(ax,'x [px]'); ylabel(ax,'y [px]');
    title(ax,'Original-image |R| with old (white) / updated (red) path');
    if panel==2, xlim(ax,bounds(:,1)'); ylim(ax,bounds(:,2)'); end
end
exportgraphics(gcf,fullfile(folder,'RethorePathComparison.png'),'Resolution',160);
% Refresh the requested final-displacement initialization and raw comparisons.
comparePhaseFieldDIC;
