clc
clear all
close all
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
addpath(correli_pathdef)
set(groot,...
'defaulttextinterpreter','latex',...
'defaultAxesTickLabelInterpreter','latex',...
'defaultLegendInterpreter','latex',...
'defaultAxesFontSize',20, ...
'defaultTextFontSize',20)


%% Import images

f = imread(fullfile(projectRoot,'Slike','reference.tif'));
g = imread(fullfile(projectRoot,'Slike','deformed.tif'));
[n,m] = size(f);

%% Show images

fig = figure('Color', 'w', 'Position', [80, 80, 1100, 500]);

% Reference image
subplot(121)
colormap gray
imagesc(f)
axis image
xlabel("$x$ [px]")
ylabel("$y$ [px]")
title("Reference")
% Deformed image
subplot(122)
colormap gray
imagesc(g)
axis image
xlabel("$x$ [px]")
ylabel("$y$ [px]")
title("Deformed")

f = single(f);
g = single(g);

%% Reading the mesh
addpath(fullfile(projectRoot,'Mesh'))
CrackMesh1;
Mesh.pos = Mesh.pos+[m/2,n/2-90,0];



frgb = repmat((f-min(f(:)))/(max(f(:)-min(f(:)))), 1, 1, 3);

fig = figure('Color', 'w', 'Position', [80, 0, 1000, 1000]);
imagesc(frgb);
mesh_plot(Mesh,"Wireframe",true,"Edgecolor","m")
axis image

%% Correlate
init = zeros(size(Mesh.pos));

fc = imgaussfilt(f,2);
gc = imgaussfilt(g,2);


%% Pyramid scheme
a = correlate(fc,gc,Mesh,init,"plot_freq",0,"max_iter",500,"conv_limit",1e-6);

[a,r,R] = correlate(f,g,Mesh,a,"plot_freq",0,"max_iter",500,"conv_limit",1e-6);

%%
fig = figure('Color', 'w', 'Position', [80, 0, 1000, 1000]);
colormap gray
imagesc(R);
axis image


%% Persist regular-DIC initialization for the phase-field workflow
measurementMask = isfinite(R); % actual kernel ROI, including mesh exclusions
imageSize = size(f);
grayScale = double(max(f(:))-min(f(:)));
if grayScale <= 0, error('example2Dscreening:Range','Reference image is constant.'); end
initialFolder = fullfile(projectRoot,'InitialValues');
if ~isfolder(initialFolder), mkdir(initialFolder); end
save(fullfile(initialFolder,'DisplacementInit.mat'), ...
    'a','Mesh','R','measurementMask','imageSize','grayScale','-v7');
residualFolder = fullfile(projectRoot,'Residual');
if ~isfolder(residualFolder), mkdir(residualFolder); end
save(fullfile(residualFolder,'R.mat'),'R','-v7');