%% Automatic residual thresholding, thinning and Euclidean crack path
% Run after screenResidualImage to reuse rho, or standalone to compute it.
% An existing rho is used unchanged. Clear rho to recompute from R.mat.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
windowSize=[23 23]; % used only if rho is not already present
numBins=256;
fixedFieldSource = 'rho'; % 'rho' or 'residual' (absolute raw residual)
fixedFieldOptions = struct('maxIterations',30,'gradientGain',0.1, ...
    'maxNormalStep',0.5,'smoothLength',3,'motionTolerance',0.02);
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
preRefinementPathXY = crackPathXY;
assert(isequal(size(rho),size(data.R),initial.imageSize),'Input image sizes must agree.');
switch fixedFieldSource
    case 'rho'
        fixedField = rho;
    case 'residual'
        fixedField = abs(double(data.R));
    otherwise
        error('extractResidualCrackPath:Field','Choose rho or residual.');
end
fixedField(~initial.measurementMask | ~isfinite(data.R)) = NaN;
fixedFieldResult = refineCrackPathFixedField( ...
    preRefinementPathXY,fixedField,initial.Mesh,fixedFieldOptions);
crackPathXY = fixedFieldResult.pathXY;
crackPathMask = fixedFieldResult.pathMask;
refinedCrackGraph = fixedFieldResult.graph;
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
title(ax,'Binary mask and fixed-field refined path');
ax=nexttile;
imagesc(ax,crackSkeleton);
axis(ax,'image');
title(ax,'One-pixel-wide skeleton');
ax=nexttile;
plot(ax,preRefinementPathXY(:,1),preRefinementPathXY(:,2),'k:');
hold(ax,'on');
plot(ax,refinedCrackGraph,'XData',refinedCrackGraph.Nodes.X, ...
    'YData',refinedCrackGraph.Nodes.Y,'Marker','none','NodeLabel',{},'EdgeColor','r');
axis(ax,'image');
set(ax,'YDir','reverse');
title(ax,'Initial (dotted) and fixed-field refined graph (red)');
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
    'preRefinementPathXY','refinedCrackGraph','fixedFieldSource','-v7');
save(fullfile(crackFolder,'FixedFieldRefinement.mat'), ...
    'fixedFieldResult','fixedFieldSource','preRefinementPathXY','-v7.3');
exportgraphics(gcf,fullfile(crackFolder,'CrackExtraction.png'),'Resolution',160);

% A larger standalone overlay makes pixel-level alignment easier to inspect.
overlayFigure = figure('Color','w','Position',[50 50 1400 400]);
ax = axes(overlayFigure);
plotCrackOnDeformed(ax,g,crackPathReferenceXY,initial.Mesh,initial.a, ...
    'Extracted crack mapped by initial displacement onto original deformed image');
exportgraphics(overlayFigure,fullfile(crackFolder,'CrackOnDeformed.png'),'Resolution',160);

% Save the actual accepted progression against the same frozen score field.
progressFigure = figure('Color','w','Position',[50 50 1200 500]);
tiledlayout(1,2);
ax = nexttile;
h = imagesc(ax,fixedField);
h.AlphaData = isfinite(fixedField);
colormap(ax,flipud(gray));
axis(ax,'image');
hold(ax,'on');
selected = unique(round(linspace(1,numel(fixedFieldResult.pathHistory),4)));
colors = lines(numel(selected));
for k = 1:numel(selected)
    points = fixedFieldResult.pathHistory{selected(k)};
    plot(ax,points(:,1),points(:,2),'Color',colors(k,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('Step %d',selected(k)-1));
end
legend(ax,'Location','best');
title(ax,'Gradient shape evolution on fixed field');
ax = nexttile;
plot(ax,0:numel(fixedFieldResult.scoreHistory)-1,fixedFieldResult.scoreHistory,'o-');
xlabel(ax,'Accepted shape step');
ylabel(ax,'Mean fixed-field score [GL]');
grid(ax,'on');
exportgraphics(progressFigure,fullfile(crackFolder,'FixedFieldRefinement.png'),'Resolution',160);
