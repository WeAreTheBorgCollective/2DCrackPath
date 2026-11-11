function result = evolveXDIC(f,g,Mesh,path,initial,folder,options)
%EVOLVEXDIC Alternate residual-gradient geometry proposals and X-DIC solves.
% Warm starts transfer enriched amplitudes by node ID after every prep change.
% Accept only lower re-solved RMS on the initial fixed pixel mask.
if nargin<7
    options = struct;
end
o = struct('iterations',8,'maxBacktracks',4,'minimumImprovement',1e-6, ...
    'levelSet',struct,'correlation',struct('maxIterations',100));
names = fieldnames(options);
for k = 1:numel(names)
    if ~isfield(o,names{k})
        error('xdic:Option','Unknown evolution option %s',names{k});
    end
    o.(names{k}) = options.(names{k});
end
validateattributes(o.iterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(o.maxBacktracks,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(o.minimumImprovement,{'numeric'},{'scalar','finite','nonnegative'});
if ~isfolder(folder)
    mkdir(folder);
end
LS = crackLevelSets2D(path,size(f));
[solution,~,R,stat] = correlateXDIC(f,g,Mesh,path,initial,o.correlation);
fixedMask = isfinite(R);
rmsValue = sqrt(mean(R(fixedMask).^2));
history = rmsValue;
paths = {path};
info = struct('accepted',true,'reason','initial');
saveXDICEvolutionStep(folder,0,LS,solution,R,stat.prep,info,g,Mesh);
status = 'maximum-iterations';
attempts = {};
for iteration = 1:o.iterations
    accepted = false;
    for backtrack = 0:o.maxBacktracks
        lsOptions = o.levelSet;
        if ~isfield(lsOptions,'duration')
            lsOptions.duration = 1;
        end
        lsOptions.duration = lsOptions.duration/2^backtrack;
        [candidate,info] = updateXDICLevelSets(LS,R,fixedMask,lsOptions);
        info.outerIteration = iteration;
        info.outerBacktrack = backtrack;
        if ~info.accepted
            attempts{end+1} = info; %#ok<AGROW>
            continue
        end
        prep = prepareXDIC(Mesh,size(f),candidate.pathXY);
        coefficients = zeros(size(prep.basis,2),2);
        nn = size(Mesh.pos,1);
        coefficients(1:nn,:) = solution.standard;
        [present,oldIndex] = ismember(prep.enrichedNodes,solution.enrichedNodes);
        coefficients(nn+find(present),:) = solution.enriched(oldIndex(present),:);
        warm = struct('coefficients',coefficients);
        [trial,~,trialR,trialStat] = correlateXDIC(f,g,Mesh, ...
            candidate.pathXY,warm,o.correlation);
        coverage = all(isfinite(trialR(fixedMask)));
        trialRMS = sqrt(mean(trialR(fixedMask).^2));
        info.trialRMS = trialRMS;
        info.previousRMS = rmsValue;
        info.solveStatus = trialStat.status;
        info.accepted = coverage && trialRMS<rmsValue-o.minimumImprovement;
        attempts{end+1} = info; %#ok<AGROW>
        if ~info.accepted
            continue
        end
        LS = candidate;
        solution = trial;
        R = trialR;
        stat = trialStat;
        rmsValue = trialRMS;
        history(end+1,1) = rmsValue; %#ok<AGROW>
        paths{end+1} = LS.pathXY; %#ok<AGROW>
        accepted = true;
        fprintf('Geometry %d: RMS %.6f, motion %.4f px, enriched nodes %d\n', ...
            iteration,rmsValue,info.maxPathMotion,numel(solution.enrichedNodes));
        saveXDICEvolutionStep(folder,iteration,LS,solution,R,stat.prep,info,g,Mesh);
        break
    end
    if ~accepted
        status = 'no-improving-geometry';
        break
    end
end
result = struct('levelSets',LS,'solution',solution,'R',R,'stat',stat, ...
    'rmsHistory',history,'paths',{paths},'status',status,'options',o, ...
    'attempts',{attempts},'fixedMask',fixedMask);
save(fullfile(folder,'Evolution.mat'),'result','-v7.3');
end
