function [solution,r,R,stat] = correlateXDIC(f,g,Mesh,crackPathXY,initial,options)
%CORRELATEXDIC MATLAB solve with standalone 2D shifted-Heaviside MEX kernels.
% The reference-gradient tangent matches Correli's modified Gauss-Newton form.
% Coefficients use [ux1 uy1 ux2 uy2 ...], standard nodes then enriched nodes.
% No BCC, mechanics, phase-field damage or crack-tip functions are included.
if nargin<6
    options = struct;
end
defaults = struct('maxIterations',100,'tolerance',1e-6, ...
    'damping',1e-8,'interpolation','cubic','maxBacktracks',15);
names = fieldnames(options);
for k = 1:numel(names)
    if ~isfield(defaults,names{k})
        error('xdic:Option','Unknown option: %s',names{k});
    end
    defaults.(names{k}) = options.(names{k});
end
p = defaults;
validateattributes(p.maxIterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(p.tolerance,{'numeric'},{'scalar','finite','positive'});
validateattributes(p.damping,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(p.maxBacktracks,{'numeric'},{'scalar','integer','nonnegative'});
p.interpolation = validatestring(p.interpolation,{'linear','cubic'});
if ~isequal(size(f),size(g))
    error('xdic:Images','Images must have equal sizes.');
end
prep = prepareXDIC(Mesh,size(f),crackPathXY);
nb = size(prep.basis,2);
nn = prep.nNodes;
coefficients = zeros(nb,2);
if nargin>=5 && ~isempty(initial)
    if isstruct(initial)
        validateattributes(initial.coefficients,{'double'},{'finite','real','size',[nb 2]});
        coefficients = initial.coefficients;
    else
        validateattributes(initial,{'numeric'},{'finite','real','2d','nrows',nn});
        if size(initial,2)<2
            error('xdic:Initial','Initialization must contain ux and uy.');
        end
        coefficients(1:nn,:) = double(initial(:,1:2));
    end
end
a = reshape(coefficients.',[],1);
state = xdicCorrMex(prep,f,g,a,p.interpolation,true);
fixedMask = state.validMask;
if nnz(fixedMask)==0
    error('xdic:Coverage','No valid initial correlation pixels.');
end
objective = mean(state.R(fixedMask).^2);
history = objective;
status = 'maximum-iterations';
for iteration = 1:p.maxIterations
    diagonal = full(diag(state.H));
    active = find(diagonal>max(diagonal)*1e-12 & diagonal>0);
    if isempty(active)
        status = 'unobservable';
        break
    end
    scale = 1./sqrt(diagonal(active));
    S = spdiags(scale,0,numel(active),numel(active));
    scaledH = S*state.H(active,active)*S;
    da = zeros(size(a));
    da(active) = scale.*((scaledH+p.damping*speye(numel(active)))\(scale.*state.b(active)));
    accepted = false;
    for backtrack = 0:p.maxBacktracks
        alpha = 2^(-backtrack);
        trial = xdicCorrMex(prep,f,g,a+alpha*da,p.interpolation,false);
        if all(trial.validMask(fixedMask))
            trialObjective = mean(trial.R(fixedMask).^2);
            if trialObjective<=objective
                accepted = true;
                break
            end
        end
    end
    if ~accepted
        status = 'line-search-stalled';
        break
    end
    a = a+alpha*da;
    objective = trialObjective;
    history(end+1,1) = objective; %#ok<AGROW>
    % Convergence is measured in the actual pixel displacement basis.
    pixelIncrement = prep.basis*reshape(alpha*da,2,[])';
    maxIncrement = max(vecnorm(pixelIncrement,2,2));
    if maxIncrement<=p.tolerance
        status = 'converged';
        break
    end
    state = xdicCorrMex(prep,f,g,a,p.interpolation,true);
end
state = xdicCorrMex(prep,f,g,a,p.interpolation,false);
coefficients = reshape(a,2,[])';
solution = struct('coefficients',coefficients,'standard',coefficients(1:nn,:), ...
    'enriched',coefficients(nn+1:end,:),'enrichedNodes',prep.enrichedNodes, ...
    'U',state.U,'coordinateFrame','reference');
R = state.R;
r = sqrt(mean(R(fixedMask).^2));
stat = struct('status',status,'iterations',iteration,'objectiveHistory',history, ...
    'fixedMask',fixedMask,'prep',prep,'options',p,'activeDOFs',active);
end
