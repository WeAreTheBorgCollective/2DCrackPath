function saveLevelSetStep(folder,step,levelSets,Damage,a,R,measurementMask,options,stage,g,Mesh)
%SAVELEVELSETSTEP Save reproducible geometry and its original-image residual.
% R belongs to the accepted a at this step, with no image damage weighting.
screenInput = double(R);
screenInput(~measurementMask) = NaN;
[screenedResidual,~,~] = screenResidual2D( ...
    screenInput,options.screeningWindowSize,0,0,0);
filename = fullfile(folder,sprintf('Step_%02d',step));
crackPathReferenceXY = levelSets.pathXY;
[crackPathDeformedXY,crackPathDisplacementXY] = mapCrackPathToDeformed( ...
    crackPathReferenceXY,Mesh,a);
coordinateFrame = 'reference';
mappingConvention = 'x = X + u(X)';
save([filename '.mat'],'step','stage','levelSets','Damage','a','R', ...
    'screenedResidual','measurementMask','options', ...
    'crackPathReferenceXY','crackPathDeformedXY','crackPathDisplacementXY', ...
    'coordinateFrame','mappingConvention','-v7');

fig = figure('Color','w','Visible','off','Position',[50 50 1300 800]);
cleanup = onCleanup(@()close(fig)); %#ok<NASGU>
tiledlayout(fig,3,2,'TileSpacing','compact');
path = levelSets.pathXY;
bounds = [min(path,[],1)-20; max(path,[],1)+20];
fields = {abs(double(R)),screenedResidual,levelSets.phi,Damage};
titles = {'Undamaged-image |f - g(X+uPF)| [GL]', ...
    sprintf('Screening %d x %d [GL]',options.screeningWindowSize), ...
    'Crack level set phi and tip zero contours','Damage'};
for panel = 1:4
    ax = nexttile;
    im = imagesc(ax,fields{panel});
    im.AlphaData = measurementMask & isfinite(fields{panel});
    axis(ax,'image');
    hold(ax,'on');
    plot(ax,path(:,1),path(:,2),'r-','LineWidth',1.5);
    if panel == 3
        contour(ax,levelSets.psiStart,[0 0],'w--');
        contour(ax,levelSets.psiEnd,[0 0],'w--');
    end
    xlim(ax,bounds(:,1)');
    ylim(ax,bounds(:,2)');
    colorbar(ax);
    title(ax,titles{panel});
    xlabel(ax,'x [px]');
    ylabel(ax,'y [px]');
end
ax = nexttile([1 2]);
plotCrackOnDeformed(ax,g,crackPathReferenceXY,Mesh,a, ...
    'Crack mapped by this step displacement onto original deformed image');
sgtitle(fig,sprintf('Step %d: %s',step,stage));
exportgraphics(fig,[filename '.png'],'Resolution',140);
end
