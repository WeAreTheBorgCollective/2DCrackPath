% Residual-gradient crack optimization after reference-frame crack extraction.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
addpath(correli_pathdef);
initial = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
data = load(fullfile(projectRoot,'CrackPaths','CrackPathInit.mat'));
if isfield(data,'crackPathReferenceXY')
    path = data.crackPathReferenceXY;
elseif isfield(data,'coordinateFrame') && strcmp(data.coordinateFrame,'reference')
    path = data.crackPathXY;
else
    error('xdic:Frame','Reference-frame crack coordinates are required.');
end
f = single(imread(fullfile(projectRoot,'Slike','reference.tif')));
g = single(imread(fullfile(projectRoot,'Slike','deformed.tif')));
options = struct;
options.iterations = 8;
% Independent of extraction screening. Set [1 1] for the raw |R| gradient.
options.levelSet = struct('screeningWindowSize',[5 5], ...
    'gradientGain',0.1,'maxNormalStep',1,'smoothLength',5);
folder = fullfile(projectRoot,'XDIC','Evolution');
result = evolveXDIC(f,g,initial.Mesh,path,initial.a,folder,options);
fprintf('Evolution status: %s; initial RMS %.6f, final %.6f\n', ...
    result.status,result.rmsHistory(1),result.rmsHistory(end));
fig = figure('Color','w','Position',[50 50 1400 700]);
tiledlayout(2,1);
ax = nexttile;
imagesc(ax,f);
colormap(ax,gray);
axis(ax,'image');
hold(ax,'on');
colors = turbo(numel(result.paths));
for k = 1:numel(result.paths)
    p = result.paths{k};
    plot(ax,p(:,1),p(:,2),'Color',colors(k,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('Step %d',k-1));
end
xlim(ax,[min(path(:,1))-15 max(path(:,1))+15]);
ylim(ax,[min(path(:,2))-15 max(path(:,2))+15]);
legend(ax,'Location','eastoutside');
title(ax,'Crack evolution in reference coordinates');
ax = nexttile;
plot(ax,0:numel(result.rmsHistory)-1,result.rmsHistory,'o-','LineWidth',1.5);
xlabel(ax,'Accepted geometry step');
ylabel(ax,'Residual RMS [GL]');
grid(ax,'on');
title(ax,'Original-image residual after each XDIC re-solve, fixed pixel mask');
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,fullfile(folder,'EvolutionSummary.png'),'Resolution',140);
