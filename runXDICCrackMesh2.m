%% Post-process the last accepted geometry from runXDICCrackEvolution.
% Historical filename retained. No new mesh, correlation or evolution is run.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
if ~exist('xdicEvolutionFolder','var')
    xdicEvolutionFolder = fullfile(projectRoot,'XDIC','Evolution');
end
sourceFile = fullfile(xdicEvolutionFolder,'Evolution.mat');
savedEvolution = load(sourceFile,'result');
result = savedEvolution.result;
[U,crackPathReferenceXY] = evaluateXDICDisplacement(result);
ux = U(:,:,1);
uy = U(:,:,2);
magnitude = hypot(ux,uy);
coordinateFrame = 'reference';
save(fullfile(xdicEvolutionFolder,'FinalDisplacement.mat'), ...
    'U','ux','uy','magnitude','crackPathReferenceXY','coordinateFrame','sourceFile','-v7.3');
fig = figure('Color','w','Position',[40 40 1500 600], ...
    'Name','Final XDIC displacement and crack geometry');
tiledlayout(fig,1,3,'TileSpacing','compact');
fields = {ux,uy,magnitude};
titles = {'Horizontal displacement [px]','Vertical displacement [px]', ...
    'Displacement magnitude [px]'};
for component = 1:3
    ax = nexttile;
    h = imagesc(ax,fields{component});
    h.AlphaData = isfinite(fields{component});
    axis(ax,'image');
    colormap(ax,turbo);
    colorbar(ax);
    hold(ax,'on');
    plot(ax,crackPathReferenceXY(:,1),crackPathReferenceXY(:,2), ...
        'k-','LineWidth',1.5);
    xlabel(ax,'Reference x [px]');
    ylabel(ax,'Reference y [px]');
    title(ax,titles{component});
end
sgtitle(fig,'Final accepted crack (black); enriched displacement on the reference grid');
exportgraphics(fig,fullfile(xdicEvolutionFolder,'FinalDisplacement.png'),'Resolution',160);
fprintf('Reconstructed displacement from %s\n',sourceFile);
