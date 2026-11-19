% Compare residual segmentation without overwriting CrackPaths/CrackPathInit.mat.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
data = load(fullfile(projectRoot,'Residual','R.mat'),'R');
initial = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
g = imread(fullfile(projectRoot,'Slike','deformed.tif'));
assert(isequal(size(data.R),size(g),initial.imageSize),'Image sizes must agree.');
raw = abs(double(data.R));
[rho,~,~] = screenResidual2D(raw,[23 23],0,0,0);
scores = {rho,raw,rho};
names = {'ValleyOtsu_screen23','Sauvola_raw','Sauvola_screen23'};
sauvolaOptions = struct('windowSize',[51 51],'k',0.2, ...
    'dynamicRange',0.5,'scalePercentile',99.5,'minCoverage',0.5);
% Short bridges only: a gap is inferred geometry, not measured crack evidence.
graphOptions = struct('MaxGap',25,'MinComponentDiameter',25,'MaxVertexDensity',2);
root = fullfile(projectRoot,'CrackPaths','Sauvola');
if ~isfolder(root)
    mkdir(root);
end
runFolder = fullfile(root,['Run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
mkdir(runFolder);
metrics = zeros(3,7);
endpoints = nan(3,4);
for trial = 1:3
    score = scores{trial};
    if trial==1
        [binaryCrack,threshold,thresholdInfo] = valleyOtsuThreshold(score,256);
    else
        [binaryCrack,threshold,thresholdInfo] = sauvolaResidualThreshold(score,sauvolaOptions);
    end
    crackSkeleton = thinCrackMask(binaryCrack);
    graphOptions.ValidMask = isfinite(score);
    [crackPathReferenceXY,crackPathMask,crackGraph,graphInfo] = ...
        connectCrackGraph(crackSkeleton,graphOptions);
    crackPathXY = crackPathReferenceXY;
    crackPathDeformedXY = zeros(0,2);
    crackPathDisplacementXY = zeros(0,2);
    inferredPathLength = 0;
    if size(crackPathXY,1)>=2
        [crackPathDeformedXY,crackPathDisplacementXY] = mapCrackPathToDeformed( ...
            crackPathXY,initial.Mesh,initial.a);
        ids = graphInfo.pathNodeIDs;
        edgeIDs = findedge(crackGraph,ids(1:end-1),ids(2:end));
        inferredPathLength = sum(crackGraph.Edges.Weight(edgeIDs).*crackGraph.Edges.IsBridge(edgeIDs));
        endpoints(trial,:) = [crackPathXY(1,:) crackPathXY(end,:)];
    end
    coordinateFrame = 'reference';
    imageSize = size(score);
    mappingConvention = 'x = X + u(X)';
    mappingDisplacementSource = 'InitialValues/DisplacementInit.mat';
    folder = fullfile(runFolder,names{trial});
    mkdir(folder);
    save(fullfile(folder,'Segmentation.mat'),'score','binaryCrack','crackSkeleton', ...
        'threshold','thresholdInfo','graphInfo','graphOptions','inferredPathLength','-v7.3');
    % An unusable trial must not look like a usable initial crack.
    if size(crackPathXY,1)>=2
        save(fullfile(folder,'CrackPathInit.mat'),'crackPathXY','crackPathReferenceXY', ...
            'crackPathDeformedXY','crackPathDisplacementXY','crackPathMask', ...
            'imageSize','coordinateFrame','mappingConvention','mappingDisplacementSource', ...
            'graphInfo','inferredPathLength','thresholdInfo','-v7');
    end
    metrics(trial,:) = [nnz(binaryCrack),graphInfo.initialComponents, ...
        graphInfo.remainingComponents,graphInfo.bridgeCount,graphInfo.pathLength, ...
        inferredPathLength,graphInfo.excludedNodeCount];
    fig = figure('Visible','off','Color','w','Position',[40 40 1400 800]);
    tiledlayout(fig,3,1,'TileSpacing','compact');
    ax = nexttile;
    im = imagesc(ax,score);
    im.AlphaData = isfinite(score);
    colormap(ax,flipud(gray(256)));
    clim(ax,[0 max(eps,max(score(isfinite(score))))]);
    axis(ax,'image');
    colorbar(ax);
    title(ax,[names{trial} ': absolute residual score [GL]'],'Interpreter','none');
    ax = nexttile;
    imagesc(ax,binaryCrack);
    colormap(ax,flipud(gray(2)));
    axis(ax,'image');
    hold(ax,'on');
    plot(ax,crackPathXY(:,1),crackPathXY(:,2),'r-','LineWidth',1.5);
    title(ax,sprintf('Binary segmentation; path %.2f px, inferred %.2f px',graphInfo.pathLength,inferredPathLength));
    ax = nexttile;
    imagesc(ax,g);
    colormap(ax,gray);
    axis(ax,'image');
    hold(ax,'on');
    plot(ax,crackPathDeformedXY(:,1),crackPathDeformedXY(:,2),'r-','LineWidth',1.5);
    title(ax,'Path mapped onto original deformed image');
    set(findall(fig,'Type','Axes'),'Toolbar',[]);
    exportgraphics(fig,fullfile(folder,'FullComparison.png'),'Resolution',140);
    % Same display-only crack-region bounds for all methods.
    reference = load(fullfile(projectRoot,'CrackPaths','CrackPathInit.mat'),'crackPathReferenceXY');
    p = reference.crackPathReferenceXY;
    for ax = findall(fig,'Type','Axes')'
        xlim(ax,[min(p(:,1))-40 max(p(:,1))+40]);
        ylim(ax,[min(p(:,2))-35 max(p(:,2))+35]);
    end
    set(findall(fig,'Type','Axes'),'Toolbar',[]);
    exportgraphics(fig,fullfile(folder,'Comparison.png'),'Resolution',140);
    close(fig);
    fprintf('%s: length %.3f, inferred %.3f, components %d, endpoints [%s]\n', ...
        names{trial},graphInfo.pathLength,inferredPathLength,graphInfo.remainingComponents,num2str(endpoints(trial,:)));
end
summary = array2table([metrics endpoints],'VariableNames', ...
    {'ForegroundPixels','InitialComponents','RemainingComponents','Bridges', ...
    'PathLength_px','InferredPathLength_px','ExcludedNodes','StartX','StartY','EndX','EndY'});
summary = addvars(summary,string(names'),'Before',1,'NewVariableNames','Method');
writetable(summary,fullfile(runFolder,'Metrics.csv'));
save(fullfile(runFolder,'Summary.mat'),'summary','sauvolaOptions','graphOptions');
save(fullfile(root,'LatestRun.mat'),'runFolder');
disp(summary);
fprintf('Trial results: %s\nExisting initial crack was not overwritten.\n',runFolder);
