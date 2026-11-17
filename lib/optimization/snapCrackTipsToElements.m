function [snapped,info] = snapCrackTipsToElements(path,Mesh,options,validMask)
% Snap an interior tip to its triangle's forward edge when remaining distance
% is strictly less than fraction times the full triangle chord along the tip
% direction. A tip already on an edge/vertex is never advanced again.
if nargin<3
    options = struct;
end
if nargin<4
    validMask = [];
end
o = struct('enabled',true,'fraction',0.25);
for name = fieldnames(options)'
    assert(isfield(o,name{1}),'Unknown tip snap option: %s',name{1});
    o.(name{1}) = options.(name{1});
end
validateattributes(o.enabled,{'logical'},{'scalar'});
validateattributes(o.fraction,{'numeric'},{'scalar','finite','positive','<=',1});
validateattributes(path,{'numeric'},{'2d','ncols',2,'finite','real'});
snapped = path;
template = struct('snapped',false,'reason','disabled','element',NaN, ...
    'remainingDistance',NaN,'directionalLength',NaN,'extension',0, ...
    'originalTip',[NaN NaN],'target',[NaN NaN]);
info = repmat(template,2,1);
if ~o.enabled || size(path,1)<2
    return
end
tri = triangulation(double(Mesh.con(:,2:4)),double(Mesh.pos(:,1:2)));
for tip = 1:2
    ids = 1:size(path,1);
    if tip==2
        ids = fliplr(ids);
    end
    ordered = path(ids,:);
    point = ordered(1,:);
    info(tip).originalTip = point;
    [element,bary] = pointLocation(tri,point);
    info(tip).element = element;
    if isnan(element)
        info(tip).reason = 'outside-mesh';
        continue
    end
    if min(bary)<=1e-9
        info(tip).reason = 'already-on-border';
        continue
    end
    next = find(vecnorm(ordered(2:end,:)-point,2,2)>1e-10,1)+1;
    if isempty(next)
        info(tip).reason = 'degenerate-path';
        continue
    end
    direction = point-ordered(next,:);
    direction = direction/norm(direction);
    xy = tri.Points(tri.ConnectivityList(element,:),:);
    B = [xy';ones(1,3)];
    slope = (B\[direction';0])';
    forward = min(-bary(slope<-1e-12)./slope(slope<-1e-12));
    backward = max(-bary(slope>1e-12)./slope(slope>1e-12));
    chord = forward-backward;
    info(tip).remainingDistance = forward;
    info(tip).directionalLength = chord;
    if isempty(forward) || isempty(backward) || forward>=o.fraction*chord-1e-10
        info(tip).reason = 'not-close-enough';
        continue
    end
    % The crack must actually enter this element, not lie wholly inside it.
    otherBary = B\[ordered(2:end,:)';ones(1,size(ordered,1)-1)];
    if ~any(min(otherBary,[],1)<=1e-9)
        info(tip).reason = 'path-wholly-inside-element';
        continue
    end
    target = point+forward*direction;
    candidate = snapped;
    % Extend the existing terminal segment collinearly; no duplicate vertex.
    if tip==1
        candidate(1,:) = target;
    else
        candidate(end,:) = target;
    end
    if ~isempty(validMask) && ~crackPathInDomain(candidate,validMask)
        info(tip).reason = 'invalid-or-outside-measurement-domain';
        continue
    end
    snapped = candidate;
    info(tip).snapped = true;
    info(tip).reason = 'snapped-to-exit-edge';
    info(tip).target = target;
    info(tip).extension = forward;
end
end
