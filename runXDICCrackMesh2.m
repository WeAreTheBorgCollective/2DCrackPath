% New-mesh example: regular DIC -> straight-crack XDIC -> level-set evolution.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
addpath(correli_pathdef);
if exist('xdicPrepMex','file')~=3
    buildXDICKernels;
end
f = single(imread(fullfile(projectRoot,'Slike','reference.tif')));
g = single(imread(fullfile(projectRoot,'Slike','deformed.tif')));
run(fullfile(projectRoot,'Mesh','CrackMesh2.m'));
meshTranslation = [size(f,2)/2,size(f,1)/2-90,0];
Mesh.pos = Mesh.pos+meshTranslation;
Mesh.MIN = min(Mesh.pos);
Mesh.MAX = max(Mesh.pos);
source = load(fullfile(projectRoot,'CrackPaths','CrackStraight.mat'));
path = double(source.crackPathReferenceXY);
tri = triangulation(Mesh.con(:,2:4),Mesh.pos(:,1:2));
assert(all(isfinite(pointLocation(tri,path))),'Straight crack must lie inside CrackMesh2.');
outputRoot = fullfile(projectRoot,'XDIC','CrackMesh2Straight');
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end
runFolder = fullfile(outputRoot,['Run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
mkdir(runFolder);
% Fresh coefficients on the new mesh; never reuse another mesh's nodal vector.
a = zeros(size(Mesh.pos));
a = correlate(imgaussfilt(f,2),imgaussfilt(g,2),Mesh,a, ...
    'plot_freq',0,'max_iter',500,'conv_limit',1e-6);
[aRegular,~,RRegular] = correlate(f,g,Mesh,a, ...
    'plot_freq',0,'max_iter',500,'conv_limit',1e-6);
save(fullfile(runFolder,'Inputs.mat'),'Mesh','meshTranslation','path', ...
    'aRegular','RRegular','-v7.3');
options = struct('iterations',12);
options.levelSet = struct('screeningWindowSize',[5 5], ...
    'gradientGain',0.1,'maxNormalStep',1,'smoothLength',5);
result = evolveXDIC(f,g,Mesh,path,aRegular,runFolder,options);
documentXDICRun(runFolder,f,g,Mesh,RRegular,result);
save(fullfile(outputRoot,'LatestRun.mat'),'runFolder');
fprintf('CrackMesh2 run: %s\nStatus: %s\n',runFolder,result.status);
