%% Automatic residual thresholding, thinning and Euclidean crack path
% Run after screenResidualImage to reuse rho, or standalone to compute it.
% An existing rho is used unchanged. Clear rho to recompute from R.mat.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
windowSize=[23 23]; % used only if rho is not already present
numBins=256;
normalSearch = struct('radius',30,'density',1); % samples per pixel on each side
tipDirectionLength = 10; % arc-length span [px] for the outward tip direction
graphOptions=struct('MaxGap',100,'MinComponentDiameter',25, ...
    'MaxVertexDensity',2);
data=load(fullfile(projectRoot,'Residual','R.mat'),'R');
initial = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
if ~exist('rho','var')
    [rho,~,~]=screenResidual2D(data.R,windowSize,0,0,0);
end
[binaryCrack,automaticThreshold,thresholdInfo]=valleyOtsuThreshold(rho,numBins);
crackSkeleton=thinCrackMask(binaryCrack);
graphOptions.ValidMask=isfinite(rho);
[crackPathXY,~,crackGraph,graphInfo]= ...
    connectCrackGraph(crackSkeleton,graphOptions);
% Extend the selected ordered skeleton path, preserving the graph for diagnostics.
skeletonPathXY = crackPathXY;
[crackPathXY,~,tipExtension] = extendCrackPathToMask( ...
    skeletonPathXY,binaryCrack,tipDirectionLength);
% Search raw |R| after endpoint extension, then replace the ordered EG nodes.
preNormalPathXY = crackPathXY;
assert(isequal(size(data.R),size(rho),initial.imageSize), ...
    'Residual, screening and initialization sizes must match.');
[crackPathXY,refinedCrackGraph,crackPathMask,normalSearchInfo] = ...
    refineCrackPathNormalMaximum(preNormalPathXY,data.R,normalSearch,initial.Mesh);
if ~normalSearchInfo.validForLevelSets
    warning('extractResidualCrackPath:NormalSearchGeometry', ...
        'Independent maxima created a crossing/collapsed/out-of-ROI path. Inspect the saved trial before XDIC.');
end
fprintf('Tip extensions: %.3f / %.3f px.\n',tipExtension.extensionLengths);
fprintf('Valley threshold: %.6g GL; foreground: %d pixels.\n', ...
    automaticThreshold,nnz(binaryCrack));
fprintf('Components: %d initial, %d removed, %d remaining; %d gap links.\n', ...
    graphInfo.initialComponents,graphInfo.removedComponents, ...
    graphInfo.remainingComponents,graphInfo.bridgeCount);
fprintf('Principal path length: %.3f px; %d graph nodes outside path.\n', ...
    graphInfo.pathLength,graphInfo.excludedNodeCount);
if graphInfo.remainingComponents>1
    warning('extractResidualCrackPath:Disconnected', ...
        'Gap limits leave %d components; returning the longest principal path.', ...
        graphInfo.remainingComponents);
end
figure('Color','w','Position',[50 50 1500 850], ...
    'Name','Automatic residual crack path');
tiledlayout(2,3,'TileSpacing','compact');
ax=nexttile;
im=imagesc(ax,rho);
im.AlphaData=isfinite(rho);
axis(ax,'image');
colorbar(ax);
title(ax,'Screened mean |R| [GL]');
ax=nexttile;
plot(ax,thresholdInfo.levels,thresholdInfo.probability);
xline(ax,automaticThreshold,'r--');
grid(ax,'on');
xlabel(ax,'Mean |R| [GL]');
ylabel(ax,'Probability');
title(ax,'Valley-emphasis threshold');
ax=nexttile;
imagesc(ax,binaryCrack);
axis(ax,'image');
hold(ax,'on');
plot(ax,crackPathXY(:,1),crackPathXY(:,2),'r-','LineWidth',1.5);
title(ax,'Binary mask and normal-refined path');
ax=nexttile;
imagesc(ax,crackSkeleton);
axis(ax,'image');
title(ax,'One-pixel-wide skeleton');
ax=nexttile;
plot(ax,refinedCrackGraph,'XData',refinedCrackGraph.Nodes.X, ...
    'YData',refinedCrackGraph.Nodes.Y,'Marker','none','NodeLabel',{}, ...
    'EdgeColor','r');
hold(ax,'on');
plot(ax,preNormalPathXY(:,1),preNormalPathXY(:,2),'k:');
axis(ax,'image');
set(ax,'YDir','reverse');
title(ax,'Normal-refined EG (red), original path (dotted)');
%% Map the extracted reference-grid path to the original deformed image
if size(crackPathXY,1)<2
    error('extractResidualCrackPath:Empty','No usable crack path; initialization was not saved.');
end
[~,g] = loadDICImages(projectRoot,initial);
if ~isequal(size(rho),size(g),initial.imageSize)
    error('extractResidualCrackPath:ImageSize','Screening, displacement and image sizes must match.');
end
crackPathReferenceXY = crackPathXY;
ax = nexttile;
[crackPathDeformedXY,crackPathDisplacementXY] = plotCrackOnDeformed( ...
    ax,g,crackPathReferenceXY,initial.Mesh,initial.a, ...
    'Mapped crack on original deformed image');

%% Persist the ordered initial crack path for level sets and phase-field DIC
if size(crackPathXY,1)<2
    error('extractResidualCrackPath:Empty','No usable crack path; initialization was not saved.');
end
crackFolder=fullfile(projectRoot,'CrackPaths');
if ~isfolder(crackFolder)
    mkdir(crackFolder);
end
imageSize=size(rho);
coordinateFrame = 'reference';
mappingConvention = 'x = X + u(X)';
mappingDisplacementSource = 'InitialValues/DisplacementInit.mat';
save(fullfile(crackFolder,'CrackPathInit.mat'), ...
    'crackPathXY','crackPathReferenceXY','crackPathDeformedXY', ...
    'crackPathDisplacementXY','coordinateFrame','mappingConvention', ...
    'mappingDisplacementSource','crackPathMask','imageSize','graphInfo', ...
    'skeletonPathXY','binaryCrack','crackSkeleton','tipExtension', ...
    'preNormalPathXY','refinedCrackGraph','crackGraph','normalSearchInfo','-v7');
exportgraphics(gcf,fullfile(crackFolder,'CrackExtraction.png'),'Resolution',160);

% A larger standalone overlay makes pixel-level alignment easier to inspect.
overlayFigure = figure('Color','w','Position',[50 50 1400 400]);
ax = axes(overlayFigure);
plotCrackOnDeformed(ax,g,crackPathReferenceXY,initial.Mesh,initial.a, ...
    'Extracted crack mapped by initial displacement onto original deformed image');
exportgraphics(overlayFigure,fullfile(crackFolder,'CrackOnDeformed.png'),'Resolution',160);
