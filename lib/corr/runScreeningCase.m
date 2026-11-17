function result = runScreeningCase(projectRoot,caseDefinition,varargin)
% Fixed-reference sequence DIC; previous total displacement seeds the next step.
p = inputParser;
addParameter(p,'ShowPlots',true,@islogical);
addParameter(p,'PublishLatest',true,@islogical);
parse(p,varargin{:});
caseName = caseDefinition.name;
imageList = caseDefinition.images;
assert(numel(imageList)>=2,'A case needs a reference and at least one deformed image.');
referenceImage = imageList{1};
f = single(imread(fullfile(projectRoot,referenceImage)));
assert(ismatrix(f),'Only grayscale 2D images are supported.');
imageSize = size(f);
grayScale = double(max(f(:))-min(f(:)));
assert(grayScale>0,'Reference image is constant.');
% Validate every input before beginning or replacing active outputs.
for k = 2:numel(imageList)
    info = imfinfo(fullfile(projectRoot,imageList{k}));
    assert(info.Height==imageSize(1) && info.Width==imageSize(2), ...
        'Image dimensions do not match: %s',imageList{k});
end
Mesh = loadScreeningMesh(fullfile(projectRoot,caseDefinition.meshFile), ...
    imageSize,caseDefinition.meshOffset);
runFolder = fullfile(projectRoot,'Results','Screening',caseName, ...
    ['Run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
mkdir(runFolder);
save(fullfile(runFolder,'Case.mat'),'caseDefinition','imageList','Mesh');
a = zeros(size(Mesh.pos));
fc = imgaussfilt(f,2);
result = struct('folder',runFolder,'completedSteps',0);
if p.Results.ShowPlots
    imageFigure = figure('Color','w','Position',[80 80 1100 500]);
    meshFigure = figure('Color','w','Position',[80 0 1000 1000]);
    frgb = repmat((f-min(f(:)))/grayScale,1,1,3);
    imagesc(frgb);
    hold on;
    mesh_plot(Mesh,'Wireframe',true,'Edgecolor','m');
    axis image;
    hold off;
    exportgraphics(meshFigure,fullfile(runFolder,'Mesh.png'));
    residualFigure = figure('Color','w','Position',[80 0 1000 1000]);
end
for stepIndex = 1:numel(imageList)-1
    deformedImage = imageList{stepIndex+1};
    g = single(imread(fullfile(projectRoot,deformedImage)));
    assert(ismatrix(g),'Only grayscale 2D images are supported.');
    initialGuess = a;
    fprintf('%s: step %d/%d, %s -> %s\n', ...
        caseName,stepIndex,numel(imageList)-1,referenceImage,deformedImage);
    gc = imgaussfilt(g,2);
    [coarse,~,~,coarseStatus] = correlate(fc,gc,Mesh,initialGuess, ...
        'plot_freq',0,'verbose_freq',-1,'max_iter',500,'conv_limit',1e-6);
    [a,r,R,status] = correlate(f,g,Mesh,coarse, ...
        'plot_freq',0,'verbose_freq',-1,'max_iter',500,'conv_limit',1e-6);
    assert(all(isfinite(a(:))) && isfinite(r), ...
        'Nonfinite solution at step %d; active output was not updated.',stepIndex);
    if status.conv_err~=0
        warning('screening:Convergence','Step %d: %s; inspect the saved status.', ...
            stepIndex,status.conv_stat);
    end
    measurementMask = isfinite(R);
    assert(any(measurementMask(:)),'No valid residual pixels.');
    stepFolder = fullfile(runFolder,sprintf('Step_%03d',stepIndex));
    mkdir(stepFolder);
    initializationFile = fullfile(stepFolder,'DisplacementInit.mat');
    residualFile = fullfile(stepFolder,'R.mat');
    save(initializationFile,'a','Mesh','R','r','measurementMask','imageSize', ...
        'grayScale','initialGuess','coarseStatus','status','caseName', ...
        'caseDefinition','referenceImage','deformedImage','stepIndex','-v7');
    save(residualFile,'R','caseName','referenceImage','deformedImage','stepIndex','-v7');
    if p.Results.PublishLatest
        initialFolder = fullfile(projectRoot,'InitialValues');
        residualFolder = fullfile(projectRoot,'Residual');
        if ~isfolder(initialFolder)
            mkdir(initialFolder);
        end
        if ~isfolder(residualFolder)
            mkdir(residualFolder);
        end
        copyfile(initializationFile,fullfile(initialFolder,'DisplacementInit.mat'));
        copyfile(residualFile,fullfile(residualFolder,'R.mat'));
    end
    result.completedSteps = stepIndex;
    save(fullfile(runFolder,'Progress.mat'),'result');
    if p.Results.ShowPlots
        figure(imageFigure);
        clf(imageFigure);
        subplot(1,2,1);
        imagesc(f);
        colormap gray;
        axis image;
        xlabel('$x$ [px]','Interpreter','latex');
        ylabel('$y$ [px]','Interpreter','latex');
        title('Reference');
        subplot(1,2,2);
        imagesc(g);
        axis image;
        xlabel('$x$ [px]','Interpreter','latex');
        ylabel('$y$ [px]','Interpreter','latex');
        title(sprintf('Deformed: step %d',stepIndex));
        exportgraphics(imageFigure,fullfile(stepFolder,'Images.png'));
        figure(residualFigure);
        clf(residualFigure);
        imagesc(R);
        colormap gray;
        axis image;
        title(sprintf('Residual: step %d',stepIndex));
        exportgraphics(residualFigure,fullfile(stepFolder,'Residual.png'));
        drawnow;
    end
end
fprintf('Saved %d steps in %s\n',result.completedSteps,runFolder);
end
