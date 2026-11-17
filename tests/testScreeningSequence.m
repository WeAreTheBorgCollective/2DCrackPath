function testScreeningSequence(runFolder)
% Validate archived results of a complete real composite sequence.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setupPaths;
previous = [];
for k = 1:11
    folder = fullfile(runFolder,sprintf('Step_%03d',k));
    saved = load(fullfile(folder,'DisplacementInit.mat'));
    residual = load(fullfile(folder,'R.mat'));
    assert(saved.stepIndex==k && residual.stepIndex==k);
    assert(strcmp(saved.referenceImage,'Slike/composite00.bmp'));
    assert(strcmp(saved.deformedImage,sprintf('Slike/composite%02d.bmp',k)));
    assert(isequaln(saved.R,residual.R));
    assert(isequal(saved.measurementMask,isfinite(saved.R)));
    assert(all(isfinite(saved.a(:))));
    if k==1
        assert(~any(saved.initialGuess(:)));
    else
        assert(isequal(saved.initialGuess,previous), ...
            'The preceding total displacement was not reused exactly.');
    end
    previous = saved.a;
    [f,g] = loadDICImages(projectRoot,saved);
    assert(isequal(g,single(imread(fullfile(projectRoot,saved.deformedImage)))));
    assert(isequal(size(f),saved.imageSize));
    assert(all(saved.Mesh.pos(:,3)==0));
    assert(all(saved.Mesh.con(:,1)==0));
    fprintf('Step %02d: %s, residual %.6g\n',k,saved.status.conv_stat,saved.r);
end
fprintf('Sequence persistence, warm starts, ROI and image mapping passed.\n');
end
