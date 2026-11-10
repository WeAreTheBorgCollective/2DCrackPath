clc
clear all
close all
addpath(correli_pathdef)
set(groot,...
'defaulttextinterpreter','latex',...
'defaultAxesTickLabelInterpreter','latex',...
'defaultLegendInterpreter','latex',...
'defaultAxesFontSize',20, ...
'defaultTextFontSize',20)


%% Import images

f = imread("Slike\reference.tif");
g = imread("Slike\deformed.tif");
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
addpath("Mesh\")
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
a = correlate(fc,gc,Mesh,init,"plot_freq",0,"max_iter",500,"conv_lim",1e-6);

[a,r,R] = correlate(f,g,Mesh,a,"plot_freq",0,"max_iter",500,"conv_lim",1e-6);

%%
fig = figure('Color', 'w', 'Position', [80, 0, 1000, 1000]);
colormap gray
imagesc(R);
axis image


% save("Residual\R.mat","R")