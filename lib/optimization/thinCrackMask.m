function skeleton = thinCrackMask(BW)
%THINCRACKMASK One-pixel-wide, topology-preserving skeleton of a binary mask.
% Uses MATLAB bwskel (Image Processing Toolbox), rather than reproducing
% the four-subiteration thinning routine used by Strini & Schiavi (2022).
% No branch pruning or component deletion is done here. False padding makes
% objects touching the image border explicit foreground/background objects.
validateattributes(BW,{'logical'},{'2d','nonempty'});
if exist('bwskel','file') ~= 2
    error('thinCrackMask:Toolbox','Image Processing Toolbox (bwskel) is required.');
end
padded = false(size(BW)+2);
padded(2:end-1,2:end-1) = BW;
padded = bwskel(padded);
skeleton = padded(2:end-1,2:end-1);
end
