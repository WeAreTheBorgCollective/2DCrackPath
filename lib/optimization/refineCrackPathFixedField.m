function result = refineCrackPathFixedField(path,field,Mesh,options)
% Gradient level-set refinement against ONE frozen scalar DIC score field.
% No image correlation, residual refresh, face growth or element-tip snapping.
if nargin<4
    options = struct;
end
o = struct('maxIterations',30,'gradientGain',0.1,'maxNormalStep',0.5, ...
    'smoothLength',3,'motionTolerance',0.02,'maxBacktracks',4);
for name = fieldnames(options)'
    assert(isfield(o,name{1}),'Unknown fixed-field option: %s',name{1});
    o.(name{1}) = options.(name{1});
end
validateattributes(o.maxIterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(o.motionTolerance,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(o.maxBacktracks,{'numeric'},{'scalar','integer','nonnegative'});
fixedField = abs(double(field));
mask = isfinite(fixedField);
assert(crackPathInDomain(path,mask),'Initial path must be valid in the fixed field ROI.');
tri = triangulation(double(Mesh.con(:,2:4)),double(Mesh.pos(:,1:2)));
assert(all(isfinite(pointLocation(tri,path))),'Initial path must lie in the mesh.');
LS = crackLevelSets2D(path,size(field));
score = scoreCrackPathInROI(path,fixedField,mask);
assert(isfinite(score),'Initial curve has unsupported score samples.');
scores = score;
paths = {path};
attempts = {};
status = 'maximum-iterations';
for iteration = 1:o.maxIterations
    accepted = false;
    for backtrack = 0:o.maxBacktracks
        % [1 1] ensures a supplied rho is NOT screened a second time.
        lsOptions = struct('velocityMethod','gradient','screeningWindowSize',[1 1], ...
            'gradientGain',o.gradientGain,'maxNormalStep',o.maxNormalStep, ...
            'smoothLength',o.smoothLength,'duration',2^(-backtrack));
        [proposal,info] = updateXDICLevelSets(LS,fixedField,mask,lsOptions);
        info.iteration = iteration;
        info.outerBacktrack = backtrack;
        if ~info.accepted
            attempts{end+1} = info; %#ok<AGROW>
            continue
        end
        [candidate,info.meshConstraint] = constrainCrackToMesh(proposal.pathXY,Mesh);
        if ~info.meshConstraint.valid || ~crackPathInDomain(candidate,mask)
            info.accepted = false;
            info.reason = 'invalid-constrained-contour';
            attempts{end+1} = info; %#ok<AGROW>
            continue
        end
        candidateScore = scoreCrackPathInROI(candidate,fixedField,mask);
        info.meanScreeningAfter = candidateScore;
        info.accepted = isfinite(candidateScore) && candidateScore>score+1e-8*max(1,score);
        if ~info.accepted
            info.reason = 'fixed-score-not-improved';
            attempts{end+1} = info; %#ok<AGROW>
            continue
        end
        motion = max([crackProjection2D(LS.pathXY,candidate); ...
            crackProjection2D(candidate,LS.pathXY)]);
        LS = crackLevelSets2D(candidate,size(field));
        score = candidateScore;
        scores(end+1,1) = score; %#ok<AGROW>
        paths{end+1} = candidate; %#ok<AGROW>
        info.maxPathMotion = motion;
        info.reason = 'accepted-fixed-field-ascent';
        attempts{end+1} = info; %#ok<AGROW>
        accepted = true;
        fprintf('Fixed-field shape %d: mean score %.6g, motion %.4f px\n', ...
            iteration,score,motion);
        break
    end
    if ~accepted
        status = 'no-improving-contour';
        break
    end
    if motion<=o.motionTolerance
        status = 'motion-tolerance';
        break
    end
end
finalPath = LS.pathXY;
n = size(finalPath,1);
G = graph((1:n-1)',(2:n)',vecnorm(diff(finalPath),2,2),n);
G.Nodes.X = finalPath(:,1);
G.Nodes.Y = finalPath(:,2);
pathMask = false(size(field));
for k = 1:n-1
    count = max(1,ceil(max(abs(finalPath(k+1,:)-finalPath(k,:)))));
    points = round(finalPath(k,:)+linspace(0,1,count+1)'.*(finalPath(k+1,:)-finalPath(k,:)));
    pathMask(sub2ind(size(field),points(:,2),points(:,1))) = true;
end
result = struct('pathXY',finalPath,'graph',G,'pathMask',pathMask, ...
    'levelSets',LS,'pathHistory',{paths},'scoreHistory',scores, ...
    'fixedField',fixedField,'mask',mask,'status',status,'options',o,'attempts',{attempts});
end
