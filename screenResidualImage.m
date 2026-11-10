%% Screen the saved 2D DIC residual (Report, Section 4.2)
% Run this script from the project root, or use run with its full path.
% Edit the settings below. The default window is 23 x 23; the threshold
% must be calibrated for the noise and gray-level scale of your 2D data.
windowSize = [23 23];       % [height width], odd pixel counts (e.g. [11 21])
residualThreshold = 5;     % mean absolute residual threshold [GL]
maxCandidates = 5;
minCandidateDistance = 5;  % Euclidean distance between centers [px]

%% Load residual and screen every valid center
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
data = load(fullfile(projectRoot, 'Residual', 'R.mat'), 'R');
if ~isfield(data, 'R')
    error('screenResidualImage:MissingR', 'Residual/R.mat must contain R.');
end
R = data.R;
[rho, candidates, eligible] = screenResidual2D(R, windowSize, ...
    residualThreshold, maxCandidates, minCandidateDistance);
disp(candidates);
fprintf('%d eligible centers; %d candidates selected.\n', ...
    nnz(eligible), height(candidates));

%% Show residual, screening map and ranked candidates
figure('Color','w','Name','2D residual screening', ...
    'Position',[80 80 1500 550]);
tiledlayout(1,3,'TileSpacing','compact');
ax = nexttile;
im = imagesc(ax,R);
im.AlphaData = isfinite(R);
axis(ax,'image'); colorbar(ax); colormap(ax,gray);
title(ax,'Residual R [GL]'); xlabel(ax,'x [px]'); ylabel(ax,'y [px]');

ax = nexttile;
im = imagesc(ax,rho);
im.AlphaData = isfinite(rho);
axis(ax,'image'); colorbar(ax);
title(ax,sprintf('Mean |R|, %d x %d px',windowSize));
xlabel(ax,'x [px]'); ylabel(ax,'y [px]');
hold(ax,'on');
for k = 1:height(candidates)
    x = candidates.x(k); y = candidates.y(k);
    rectangle(ax,'Position',[x-windowSize(2)/2, y-windowSize(1)/2, ...
        windowSize(2), windowSize(1)],'EdgeColor','w');
    plot(ax,x,y,'r+','MarkerSize',9,'LineWidth',1.5);
    text(ax,x+1,y,sprintf('%d',k),'Color','w','FontWeight','bold');
end

ax = nexttile;
bar(ax,candidates.rank,candidates.meanAbsoluteResidual);
hold(ax,'on');
yline(ax,residualThreshold,'r--','Threshold');
xlabel(ax,'Candidate rank'); ylabel(ax,'Mean |R| [GL]');
title(ax,'Ranked candidates'); grid(ax,'on');
pbaspect(ax,[1 1 1]);
% rho, eligible and candidates remain in the workspace for later steps.
