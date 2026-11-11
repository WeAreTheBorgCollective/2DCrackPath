function savePaperPFIteration(folder,iteration,model,U,d,r,Phi,history,p,gOriginal,plotBounds)
%SAVEPAPERPFITERATION Save independent damage and display in deformed x frame.
Damage = zeros(model.imageSize);
Damage(model.ids) = d;
residual = nan(model.imageSize);
residual(model.ids) = r;
drivingForce = nan(model.imageSize);
drivingForce(model.ids) = Phi;
pixelDisplacement = model.N*U;
measurementMask = model.mask;
coordinateFrame = 'deformed';
displacementConvention = 'X = x - u(x)';
damageMask = Damage>=p.damageThreshold & measurementMask;
crackSkeleton = thinCrackMask(damageMask);
graphOptions = struct('MaxGap',0,'MinComponentDiameter',5, ...
    'MaxVertexDensity',Inf,'ValidMask',measurementMask);
[crackPathDeformedXY,~,~,crackInfo] = connectCrackGraph(crackSkeleton,graphOptions);
nodes = model.nodes;
filename = fullfile(folder,sprintf('Iteration_%02d',iteration));
stateKind = 'staggered iteration';
if isempty(history) && iteration>0
    stateKind = 'fixed-displacement calibration trial';
end
save([filename '.mat'],'iteration','Damage','residual','drivingForce', ...
    'pixelDisplacement','U','nodes','measurementMask','coordinateFrame', ...
    'displacementConvention','damageMask','crackSkeleton','crackPathDeformedXY', ...
    'crackInfo','history','p','stateKind','-v7');
fig = figure('Visible','off','Color','w','Position',[40 40 1500 600]);
cleanup = onCleanup(@()close(fig)); %#ok<NASGU>
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for panel = 1:3
    if panel==3
        ax = nexttile([1 2]);
    else
        ax = nexttile;
    end
    if panel==2
        im = imagesc(ax,Damage,[0 1]);
        im.AlphaData = measurementMask;
        colormap(ax,turbo(256));
        colorbar(ax);
        title(ax,'Independent damage PDE solution d(x)');
    else
        imagesc(ax,gOriginal);
        colormap(ax,gray(256));
        grayLimits = [min(gOriginal(:)) max(gOriginal(:))];
        if panel==1
            title(ax,'Original deformed image (gray levels)');
        else
            hold(ax,'on');
            if any(damageMask(:))
                contour(ax,Damage,[p.damageThreshold p.damageThreshold], ...
                    'Color',[0 0.9 1],'LineWidth',1);
                [yy,xx] = find(crackSkeleton);
                plot(ax,xx,yy,'r.','MarkerSize',3);
            end
            if ~any(damageMask(:))
                text(ax,0.02,0.95,'No damage reaches the contour threshold', ...
                    'Units','normalized','Color','r','VerticalAlignment','top');
            end
            title(ax,sprintf('Damage contour d=%.2f (cyan), skeleton (red); no coordinate warp', ...
                p.damageThreshold));
        end
    end
    if panel~=2
        clim(ax,grayLimits);
    end
    axis(ax,'image');
    xlim(ax,plotBounds(1:2));
    ylim(ax,plotBounds(3:4));
    xlabel(ax,'Deformed x [px]');
    ylabel(ax,'Deformed y [px]');
end
sgtitle(fig,sprintf('Paper PF-DIC, %s %d',stateKind,iteration));
exportgraphics(fig,[filename '.png'],'Resolution',140);
end
