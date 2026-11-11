function metrics = documentXDICRun(folder,f,g,Mesh,RRegular,result)
%DOCUMENTXDICRUN Export fixed-mask GL residual metrics and every accepted path.
count = numel(result.paths);
mask = result.fixedMask & isfinite(RRegular);
stage = ["Regular DIC";compose("XDIC step %d",(0:count-1)')];
rmse = zeros(count+1,1);
mae = rmse;
maximum = rmse;
lengthPixels = nan(count+1,1);
nEnriched = zeros(count+1,1);
residuals = cell(count+1,1);
residuals{1} = double(RRegular);
for k = 1:count
    snapshot = load(fullfile(folder,sprintf('Step_%02d.mat',k-1)));
    residuals{k+1} = snapshot.R;
    assert(all(isfinite(snapshot.R(mask))),'Comparison mask lost coverage.');
    path = snapshot.LS.pathXY;
    lengthPixels(k+1) = sum(vecnorm(diff(path),2,2));
    nEnriched(k+1) = numel(snapshot.solution.enrichedNodes);
    points = table(path(:,1),path(:,2),snapshot.minusFace(:,1),snapshot.minusFace(:,2), ...
        snapshot.plusFace(:,1),snapshot.plusFace(:,2), ...
        'VariableNames',{'ReferenceX','ReferenceY','DeformedMinusX','DeformedMinusY','DeformedPlusX','DeformedPlusY'});
    writetable(points,fullfile(folder,sprintf('CrackPath_%02d.csv',k-1)));
end
for k = 1:count+1
    values = abs(residuals{k}(mask));
    rmse(k) = sqrt(mean(values.^2));
    mae(k) = mean(values);
    maximum(k) = max(values);
end
metrics = table(stage,rmse,mae,maximum,lengthPixels,nEnriched, ...
    repmat(nnz(mask),count+1,1),'VariableNames', ...
    {'Stage','RMS_GL','MeanAbsolute_GL','MaximumAbsolute_GL','CrackLength_px','EnrichedNodes','PixelCount'});
writetable(metrics,fullfile(folder,'ResidualMetrics.csv'));
disp(metrics);
% Identical nonnegative limits; no normalization, logarithms or percentile clipping.
limits = [0 max(maximum)];
fig = figure('Visible','off','Color','w','Position',[30 30 1500 450]);
tiledlayout(fig,1,3,'TileSpacing','compact');
selected = [1 2 count+1];
for j = 1:3
    k = selected(j);
    ax = nexttile;
    h = imagesc(ax,abs(residuals{k}));
    h.AlphaData = mask;
    colormap(ax,flipud(gray(256)));
    clim(ax,limits);
    axis(ax,'image');
    xlim(ax,[min(Mesh.pos(:,1))-5 max(Mesh.pos(:,1))+5]);
    ylim(ax,[min(Mesh.pos(:,2))-5 max(Mesh.pos(:,2))+5]);
    colorbar(ax);
    title(ax,sprintf('%s: RMS %.3f GL',stage(k),rmse(k)));
    xlabel(ax,'Reference x [px]');
    ylabel(ax,'Reference y [px]');
end
sgtitle(fig,'Absolute residual [GL]; white = zero, black = maximum (shared scale)');
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,fullfile(folder,'ResidualComparison.png'),'Resolution',140);
close(fig);
fig = figure('Visible','off','Color','w','Position',[30 30 1400 850]);
tiledlayout(fig,3,1,'TileSpacing','compact');
paths = vertcat(result.paths{:});
for j = 1:3
    k = selected(j);
    ax = nexttile;
    h = imagesc(ax,abs(residuals{k}));
    h.AlphaData = mask;
    colormap(ax,flipud(gray(256)));
    clim(ax,limits);
    axis(ax,'image');
    xlim(ax,[min(paths(:,1))-10 max(paths(:,1))+10]);
    ylim(ax,[min(paths(:,2))-15 max(paths(:,2))+15]);
    colorbar(ax);
    title(ax,sprintf('%s: RMS %.3f GL',stage(k),rmse(k)));
    xlabel(ax,'Reference x [px]');
    ylabel(ax,'y [px]');
end
sgtitle(fig,'Absolute residual [GL]; white = zero, black = high (shared scale)');
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,fullfile(folder,'ResidualCrackZoom.png'),'Resolution',140);
close(fig);
fig = figure('Visible','off','Color','w','Position',[30 30 1500 900]);
tiledlayout(fig,3,1,'TileSpacing','compact');
ax = nexttile;
imagesc(ax,f);
colormap(ax,gray);
axis(ax,'image');
hold(ax,'on');
colors = turbo(count);
for k = 1:count
    path = result.paths{k};
    plot(ax,path(:,1),path(:,2),'Color',colors(k,:),'LineWidth',1.5);
end
bounds = [min(vertcat(result.paths{:}),[],1)-15;max(vertcat(result.paths{:}),[],1)+15];
xlim(ax,bounds(:,1)');
ylim(ax,bounds(:,2)');
title(ax,sprintf('Reference crack evolution: blue start to red final; %d accepted updates',count-1));
ax = nexttile;
imagesc(ax,g);
colormap(ax,gray);
axis(ax,'image');
hold(ax,'on');
for k = 1:count
    snapshot = load(fullfile(folder,sprintf('Step_%02d.mat',k-1)),'minusFace','plusFace');
    center = (snapshot.minusFace+snapshot.plusFace)/2;
    plot(ax,center(:,1),center(:,2),'Color',colors(k,:),'LineWidth',1.2);
end
plot(ax,snapshot.minusFace(:,1),snapshot.minusFace(:,2),'c-','LineWidth',1.5);
plot(ax,snapshot.plusFace(:,1),snapshot.plusFace(:,2),'m-','LineWidth',1.5);
xlim(ax,bounds(:,1)');
ylim(ax,[bounds(1,2) bounds(2,2)+15]);
title(ax,'Deformed image: evolving midpoint; final minus/plus traces cyan/magenta');
ax = nexttile;
plot(ax,0:count-1,rmse(2:end),'o-');
grid(ax,'on');
xlabel(ax,'Accepted XDIC geometry step');
ylabel(ax,'Residual RMS [GL]');
title(ax,result.status,'Interpreter','none');
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,fullfile(folder,'CrackEvolution.png'),'Resolution',140);
close(fig);
fid = fopen(fullfile(folder,'Results.md'),'w');
cleanupFile = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# CrackMesh2 / CrackStraight run\n\n');
fprintf(fid,'Mesh: %d nodes, %d triangles. Fixed comparison mask: %d pixels.\n\n',size(Mesh.pos,1),size(Mesh.con,1),nnz(mask));
fprintf(fid,'Status: %s. Accepted updates: %d.\n\n',result.status,count-1);
fprintf(fid,'| Stage | RMS [GL] | Mean absolute [GL] | Max absolute [GL] |\n|---|---:|---:|---:|\n');
for k = 1:height(metrics)
    fprintf(fid,'| %s | %.6f | %.6f | %.6f |\n',stage(k),rmse(k),mae(k),maximum(k));
end
fprintf(fid,'\nResidualComparison.png uses absolute original-image residuals with shared limits [0, %.6f] GL; high values are black. No clipping or per-panel normalization.\n',limits(2));
fprintf(fid,'\nCrackPath_XX.csv contains reference coordinates and both deformed-image traces. Step_XX.mat saves the three level sets, coefficients and signed residual; Evolution.mat also saves rejected trials. CrackEvolution.png shows paths and RMS.\n');
fprintf(fid,'\nThe mesh is translated by [imageWidth/2,imageHeight/2-90] as in example2Dscreening. The supplied crack is already in reference image coordinates and is not translated. Regular DIC uses sigma 2 then original images. XDIC evolution uses original images with a 5x5 residual screening window.\n');
end
