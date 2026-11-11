function p = paperPFParameters
%PAPERPFPARAMETERS Paper Section 2.7 and experimental Table 3 starting values.
p = struct;
p.lc = 5;
p.mu = 0.1;
p.w1 = 0.05;
p.a1 = 3e-5;
p.c1 = 1e-5;
p.w2 = 0.04;
p.a2 = 0.85;
p.c2 = 0.05;
p.residualBlock = 4;
p.grayBlock = 8;
p.displacementSpacing = 16;
p.tolerance = 0.05;
p.maxIterations = 40;
p.damageTolerance = 1e-7;
p.damageMaxIterations = 100;
p.damping = 1e-8;
p.maxDisplacementStep = 2;
p.damageThreshold = 0.5;
% Registration safeguard beyond the paper's illustrative one-step GN.
p.displacementSmoothing = [4 2 1 0];
p.innerMaxIterations = 12;
p.innerTolerance = 1e-3;
p.invertImages = true; % Paper assumes bright crack/background; our cracks are dark.
end
