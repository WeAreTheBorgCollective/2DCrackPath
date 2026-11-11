% Compare Correli, ordinary T3 MEX DIC and shifted-Heaviside T3 X-DIC.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
addpath(correli_pathdef);
if exist('xdicPrepMex','file')~=3
    buildXDICKernels;
end
initial = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
pathData = load(fullfile(projectRoot,'CrackPaths','CrackPathInit.mat'));
if isfield(pathData,'crackPathReferenceXY')
    crack = pathData.crackPathReferenceXY;
elseif isfield(pathData,'coordinateFrame') && strcmp(pathData.coordinateFrame,'reference')
    crack = pathData.crackPathXY;
else
    error('xdic:Frame','CrackPathInit must explicitly contain reference-frame coordinates.');
end
f = single(imread(fullfile(projectRoot,'Slike','reference.tif')));
g = single(imread(fullfile(projectRoot,'Slike','deformed.tif')));
Mesh = initial.Mesh;
[aRegular,~,RRegular] = correlate(f,g,Mesh,initial.a, ...
    'plot_freq',0,'max_iter',100,'conv_limit',1e-6);
options = struct('maxIterations',100,'tolerance',1e-6);
[ordinary,~,RFE,statFE] = correlateXDIC(f,g,Mesh,[],initial.a,options);
[extended,~,RXF,statXF] = correlateXDIC(f,g,Mesh,crack,initial.a,options);
% Evaluate native coefficients with the same pixel basis for displacement plots.
p = statFE.prep;
nativeState = xdicCorrMex(p,f,g,reshape(aRegular(:,1:2)',[],1),'cubic',false);
commonMask = isfinite(RRegular) & isfinite(RFE) & isfinite(RXF);
residuals = {double(RRegular),RFE,RXF};
solutions = {nativeState.U,ordinary.U,extended.U};
labels = {'Correli','T3 FE','T3 XFEM'};
rmse = cellfun(@(R) sqrt(mean(R(commonMask).^2)),residuals);
metrics = table(string(labels)',rmse','VariableNames',{'Method','ResidualRMSE'});
disp(metrics);
fprintf('Enriched nodes: %d; FE status: %s; XFEM status: %s\n', ...
    numel(extended.enrichedNodes),statFE.status,statXF.status);
out = fullfile(projectRoot,'XDIC');
if ~isfolder(out)
    mkdir(out);
end
save(fullfile(out,'Comparison.mat'),'aRegular','ordinary','extended', ...
    'RRegular','RFE','RXF','statFE','statXF','commonMask','metrics','crack','Mesh','-v7.3');
fig = figure('Color','w','Position',[50 100 1500 550]);
tiledlayout(1,3);
limit = max(cellfun(@(R) prctile(abs(R(commonMask)),99),residuals));
for k = 1:3
    nexttile;
    h = imagesc(residuals{k},[-limit limit]);
    set(h,'AlphaData',isfinite(residuals{k}));
    axis image;
    colorbar;
    title(sprintf('%s: RMSE %.3f',labels{k},rmse(k)));
    xlabel('x [px]');
    ylabel('y [px]');
end
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,fullfile(out,'ResidualComparison.png'));
fig = figure('Color','w','Position',[50 50 1500 900]);
tiledlayout(2,3);
for component = 1:2
    limits = [inf -inf];
    for k = 1:3
        v = solutions{k}(:,:,component);
        limits = [min(limits(1),min(v(commonMask))),max(limits(2),max(v(commonMask)))];
    end
    for k = 1:3
        nexttile;
        h = imagesc(solutions{k}(:,:,component),limits);
        set(h,'AlphaData',isfinite(solutions{k}(:,:,component)));
        axis image;
        colorbar;
        title(sprintf('%s: u_%d [px]',labels{k},component));
    end
end
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,fullfile(out,'DisplacementComparison.png'));
fig = figure('Color','w');
imagesc(f);
colormap gray;
axis image;
hold on;
triplot(Mesh.con(:,2:4),Mesh.pos(:,1),Mesh.pos(:,2),'Color',[.4 .7 .8]);
plot(crack(:,1),crack(:,2),'r-','LineWidth',2);
e = statXF.prep.enrichedNodes;
t = statXF.prep.tipSupportNodes;
plot(Mesh.pos(e,1),Mesh.pos(e,2),'go','MarkerFaceColor','g');
plot(Mesh.pos(t,1),Mesh.pos(t,2),'ms');
title('Reference mesh: enriched nodes (green), tip-element nodes (magenta; may also be enriched)');
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,fullfile(out,'Enrichment.png'));
