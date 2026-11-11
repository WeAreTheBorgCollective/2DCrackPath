function [d,info] = paperPFDamageSolve(K,mass,H,d,p,colors)
%PAPERPFDAMAGESOLVE Eq.33 with the physical box 0<=d<=1.
% Raw tanh Phi may be negative. Do NOT silently clip Phi or H: minimize the
% resulting box-constrained quadratic and report projected KKT stationarity.
% Bound handling and the active-set method are numerical implementation choices.
A = p.lc*K+spdiags(mass.*(p.mu/p.lc+H),0,numel(d),numel(d));
b = mass.*H;
diagonal = full(diag(A));
if any(diagonal<=0)
    error('paperPF:DamageDiagonal','Nonpositive damage diagonal; check parameters.');
end
offDiagonal = A-spdiags(diagonal,0,numel(d),numel(d));
for color = 0:1
    ids = find(colors==color);
    if nnz(offDiagonal(ids,ids))~=0
        error('paperPF:DamageColors','Damage relaxation requires independent color sets.');
    end
end
energy = @(z) 0.5*(z'*(A*z))-b'*z;
d = min(1,max(0,d));
converged = false;
for iteration = 1:p.damageMaxIterations
    % Exact bound-constrained coordinate minimization on independent sets.
    % These sweeps reduce the quadratic even when its full Hessian is indefinite.
    for sweep = 1:10
        for color = 0:1
            ids = colors==color;
            rhs = b-offDiagonal*d;
            d(ids) = min(1,max(0,rhs(ids)./diagonal(ids)));
        end
    end
    gradient = A*d-b;
    projected = d-min(1,max(0,d-gradient));
    kkt = norm(projected,Inf);
    if kkt<=p.damageTolerance
        converged = true;
        break
    end
    free = (d>1e-12 & d<1-1e-12) | (d<=1e-12 & gradient<0) | ...
        (d>=1-1e-12 & gradient>0);
    step = zeros(size(d));
    [C,flag] = chol(A(free,free));
    if flag==0
        step(free) = -(C\(C'\gradient(free)));
    else
        step = -gradient./diagonal;
    end
    oldEnergy = energy(d);
    accepted = false;
    for backtrack = 0:20
        candidate = min(1,max(0,d+2^(-backtrack)*step));
        change = candidate-d;
        if gradient'*change<0 && energy(candidate)<=oldEnergy+1e-4*(gradient'*change)
            accepted = true;
            break
        end
    end
    if ~accepted
        % Scaled projected gradient is a descent fallback for indefinite free Hessians.
        step = -gradient./diagonal;
        for backtrack = 0:25
            candidate = min(1,max(0,d+2^(-backtrack)*step));
            change = candidate-d;
            if gradient'*change<0 && energy(candidate)<=oldEnergy+1e-4*(gradient'*change)
                accepted = true;
                break
            end
        end
    end
    if ~accepted
        break
    end
    d = candidate;
end
gradient = A*d-b;
kkt = norm(d-min(1,max(0,d-gradient)),Inf);
converged = kkt<=p.damageTolerance;
info = struct('converged',converged,'iterations',iteration, ...
    'projectedKKT',kkt,'quadraticEnergy',energy(d));
end
