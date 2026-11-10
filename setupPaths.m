function setupPaths()
%SETUPPATHS Add the project's library folders to the MATLAB path.
%   Run setupPaths from the project root before using library functions.
%   Paths are resolved relative to this file, not the current folder.
%   Changes apply to the current MATLAB session.

projectRoot = fileparts(mfilename('fullpath'));
subfolders = {'mesh', 'corr', 'plot', 'optimization'};

for k = 1:numel(subfolders)
    addpath(fullfile(projectRoot, 'lib', subfolders{k}));
end
end
