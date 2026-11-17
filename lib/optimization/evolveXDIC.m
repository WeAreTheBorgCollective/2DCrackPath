function result = evolveXDIC(f,g,Mesh,path,initial,folder,options)
%EVOLVEXDIC Alternate residual-gradient geometry proposals and X-DIC solves.
% Warm starts transfer enriched amplitudes by node ID after every prep change.
% Accept only lower re-solved RMS on the initial fixed pixel mask.
if nargin<7
    options = struct;
end
o = struct('iterations',8,'maxBacktracks',4,'minimumImprovement',1e-6, ...
    'tipSnap',struct('enabled',true,'fraction',0.25), ...
    'faceGrowth',struct('enabled',true,'maxExtension',10,'tangentLength',10), ...
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
geometryPrep = prepareXDIC(Mesh,size(f),[]);
geometryMask = false(size(f));
geometryMask(geometryPrep.pixelIndex) = true;
[path,initialSnap] = snapCrackTipsToElements(path,Mesh,o.tipSnap,geometryMask);
LS = crackLevelSets2D(path,size(f));
if isstruct(initial) && isfield(initial,'standard') && isfield(initial,'enrichedNodes')
    [solution,R,stat] = solveGeometry(f,g,Mesh,path,initial,o.correlation);
else
    assert(~isstruct(initial) || ~any([initialSnap.snapped]), ...
        'Snapping changes enrichment: provide nodal initialization or a full XDIC solution.');
    [solution,~,R,stat] = correlateXDIC(f,g,Mesh,path,initial,o.correlation);
end
fixedMask = isfinite(R);
rmsValue = sqrt(mean(R(fixedMask).^2));
history = rmsValue;
paths = {path};
info = struct('accepted',true,'reason','initial','tipSnap',initialSnap);
saveXDICEvolutionStep(folder,0,LS,solution,R,stat.prep,info,g,Mesh);
status = 'maximum-iterations';
attempts = {};
acceptedStep = 0;
faceOptions = o.faceGrowth;
faceEnabled = true;
if isfield(faceOptions,'enabled')
    faceEnabled = faceOptions.enabled;
    validateattributes(faceEnabled,{'logical'},{'scalar'});
    faceOptions = rmfield(faceOptions,'enabled');
end
for iteration = 1:o.iterations
    accepted = false;
    for backtrack = 0:o.maxBacktracks
        lsOptions = o.levelSet;
        if ~isfield(lsOptions,'duration')
            lsOptions.duration = 1;
        end
        lsOptions.duration = lsOptions.duration/2^backtrack;
        [candidate,info] = updateXDICLevelSets(LS,R,fixedMask,lsOptions);
        info.stage = 'shape';
        info.outerIteration = iteration;
        info.outerBacktrack = backtrack;
        if ~info.accepted
            attempts{end+1} = info; %#ok<AGROW>
            continue
        end
        [snappedPath,info.tipSnap] = snapCrackTipsToElements( ...
            candidate.pathXY,Mesh,o.tipSnap,fixedMask);
        candidate = crackLevelSets2D(snappedPath,size(f));
        [trial,trialR,trialStat] = solveGeometry(f,g,Mesh,candidate.pathXY,solution,o.correlation);
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
        acceptedStep = acceptedStep+1;
        saveXDICEvolutionStep(folder,acceptedStep,LS,solution,R,stat.prep,info,g,Mesh);
        break
    end
    % After shape optimization, test independent length growth. Even a stalled
    % shape proposal may have converging face rays, so do not stop beforehand.
    if faceEnabled
        [grownPath,tipInfo] = proposeXDICFaceGrowth(LS,Mesh,solution,stat.prep,faceOptions);
        growthInfo = struct('stage','face-growth','outerIteration',iteration, ...
            'accepted',false,'tips',tipInfo,'reason','no-forward-intersection');
        if any([tipInfo.accepted])
            [grownPath,growthInfo.tipSnap] = snapCrackTipsToElements( ...
                grownPath,Mesh,o.tipSnap,fixedMask);
            if crackPathInDomain(grownPath,fixedMask)
                [trial,trialR,trialStat] = solveGeometry(f,g,Mesh,grownPath,solution,o.correlation);
                trialRMS = sqrt(mean(trialR(fixedMask).^2));
                growthInfo.trialRMS = trialRMS;
                growthInfo.previousRMS = rmsValue;
                growthInfo.solveStatus = trialStat.status;
                growthInfo.accepted = all(isfinite(trialR(fixedMask))) && ...
                    trialRMS<rmsValue-o.minimumImprovement;
                growthInfo.reason = 'no-residual-improvement';
                if growthInfo.accepted
                    LS = crackLevelSets2D(grownPath,size(f));
                    solution = trial;
                    R = trialR;
                    stat = trialStat;
                    rmsValue = trialRMS;
                    history(end+1,1) = rmsValue; %#ok<AGROW>
                    paths{end+1} = LS.pathXY; %#ok<AGROW>
                    accepted = true;
                    growthInfo.reason = 'accepted';
                    acceptedStep = acceptedStep+1;
                    saveXDICEvolutionStep(folder,acceptedStep,LS,solution,R, ...
                        stat.prep,growthInfo,g,Mesh);
                    fprintf('Face growth %d: RMS %.6f, tip increments %.3f / %.3f px\n', ...
                        iteration,rmsValue,tipInfo.extension);
                end
            else
                growthInfo.reason = 'extension-outside-measurement-domain';
            end
        end
        attempts{end+1} = growthInfo; %#ok<AGROW>
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

function [trial,R,stat] = solveGeometry(f,g,Mesh,path,solution,options)
prep = prepareXDIC(Mesh,size(f),path);
nn = size(Mesh.pos,1);
coefficients = zeros(size(prep.basis,2),2);
coefficients(1:nn,:) = solution.standard;
[present,oldIndex] = ismember(prep.enrichedNodes,solution.enrichedNodes);
coefficients(nn+find(present),:) = solution.enriched(oldIndex(present),:);
warm = struct('coefficients',coefficients);
[trial,~,R,stat] = correlateXDIC(f,g,Mesh,path,warm,options);
end
