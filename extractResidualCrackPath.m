%% Automatic residual thresholding, thinning and Euclidean crack path
% Run after screenResidualImage to reuse rho, or standalone to compute it.
% An existing rho is used unchanged. Clear rho to recompute from R.mat.
projectRoot=fileparts(mfilename('fullpath'));
addpath(projectRoot); setupPaths;
windowSize=[23 23]; % used only if rho is not already present
numBins=256;
graphOptions=struct('MaxGap',100,'MinComponentDiameter',25, ...
    'MaxVertexDensity',2);
if ~exist('rho','var')
    data=load(fullfile(projectRoot,'Residual','R.mat'),'R');
    [rho,~,~]=screenResidual2D(data.R,windowSize,0,0,0);
end
[binaryCrack,automaticThreshold,thresholdInfo]=valleyOtsuThreshold(rho,numBins);
crackSkeleton=thinCrackMask(binaryCrack);
graphOptions.ValidMask=isfinite(rho);
[crackPathXY,crackPathMask,crackGraph,graphInfo]= ...
    connectCrackGraph(crackSkeleton,graphOptions);
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
ax=nexttile; im=imagesc(ax,rho); im.AlphaData=isfinite(rho);
axis(ax,'image'); colorbar(ax); title(ax,'Screened mean |R| [GL]');
ax=nexttile; plot(ax,thresholdInfo.levels,thresholdInfo.probability);
xline(ax,automaticThreshold,'r--'); grid(ax,'on');
xlabel(ax,'Mean |R| [GL]'); ylabel(ax,'Probability'); title(ax,'Valley-emphasis threshold');
ax=nexttile; imagesc(ax,binaryCrack); axis(ax,'image'); title(ax,'Binary crack mask');
ax=nexttile; imagesc(ax,crackSkeleton); axis(ax,'image'); title(ax,'One-pixel-wide skeleton');
ax=nexttile;
if numnodes(crackGraph)>0
    p=plot(ax,crackGraph,'XData',crackGraph.Nodes.X,'YData',crackGraph.Nodes.Y, ...
        'Marker','none','NodeLabel',{},'EdgeColor',[0.2 0.5 0.8]);
    bridgeIDs=find(crackGraph.Edges.IsBridge);
    if ~isempty(bridgeIDs), highlight(p,'Edges',bridgeIDs,'EdgeColor','r','LineWidth',2); end
end
axis(ax,'image'); set(ax,'YDir','reverse'); title(ax,'Retained graph (red = inferred gaps)');
ax=nexttile; im=imagesc(ax,rho); im.AlphaData=isfinite(rho);
axis(ax,'image'); hold(ax,'on');
if ~isempty(crackPathXY)
    plot(ax,crackPathXY(:,1),crackPathXY(:,2),'r-','LineWidth',1.5);
end
title(ax,sprintf('Principal path: %.1f px',graphInfo.pathLength));
%% Persist the ordered initial crack path for level sets and phase-field DIC
if size(crackPathXY,1)<2
    error('extractResidualCrackPath:Empty','No usable crack path; initialization was not saved.');
end
crackFolder=fullfile(projectRoot,'CrackPaths');
if ~isfolder(crackFolder), mkdir(crackFolder); end
imageSize=size(rho);
save(fullfile(crackFolder,'CrackPathInit.mat'), ...
    'crackPathXY','crackPathMask','imageSize','graphInfo','-v7');
