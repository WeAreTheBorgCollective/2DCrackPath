function buildXDICKernels
%BUILDXDICKERNELS Compile standalone 2D X-DIC MEX kernels.
% Requires a MATLAB-supported C++ compiler (mex -setup C++).
root = fileparts(mfilename('fullpath'));
output = fullfile(root,'lib','corr');
source = fullfile(output,'xdic_src');
for name = {'xdicPrepMex','xdicCorrMex'}
    clear(name{1});
    mex('-R2018a','-DNOMINMAX', ...
        fullfile(source,[name{1} '.cpp']),'-outdir',output);
end
addpath(output);
fprintf('Built 2D X-DIC kernels in %s\n',output);
end
