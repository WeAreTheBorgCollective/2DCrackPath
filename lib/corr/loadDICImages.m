function [f,g] = loadDICImages(projectRoot,initial)
% Read the image pair associated with the active regular-DIC initialization.
% Legacy initialization files without image metadata use the original pair.
if nargin<2
    initial = load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'), ...
        'referenceImage','deformedImage');
end
referenceImage = 'Slike/reference.tif';
deformedImage = 'Slike/deformed.tif';
if isfield(initial,'referenceImage') && isfield(initial,'deformedImage')
    referenceImage = initial.referenceImage;
    deformedImage = initial.deformedImage;
end
f = single(imread(fullfile(projectRoot,referenceImage)));
g = single(imread(fullfile(projectRoot,deformedImage)));
assert(ismatrix(f) && isequal(size(f),size(g)), ...
    'Expected matching grayscale 2D images.');
end
