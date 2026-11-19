clc
clear all
close all
%% Select a case from screeningCases.m.
% Set screeningCaseName before running to select another case.
% Example: screeningCaseName = 'composite050'; example2Dscreening
if ~exist('screeningCaseName','var')
    screeningCaseName = 'composite100';
end
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
setupPaths;
addpath(correli_pathdef);
cases = screeningCases;
assert(isfield(cases,screeningCaseName),'Unknown screening case: %s',screeningCaseName);
set(groot,'defaulttextinterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesFontSize',20,'defaultTextFontSize',20);
screeningRun = runScreeningCase(projectRoot,cases.(screeningCaseName));
% Retain the familiar workspace variables for interactive inspection.
latest = load(fullfile(screeningRun.folder, ...
    sprintf('Step_%03d',screeningRun.completedSteps),'DisplacementInit.mat'));
a = latest.a;
Mesh = latest.Mesh;
R = latest.R;
r = latest.r;
[f,g] = loadDICImages(projectRoot,latest);
