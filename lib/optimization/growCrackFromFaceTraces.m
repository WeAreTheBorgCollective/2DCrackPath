function [grown,info] = growCrackFromFaceTraces(path,minusFace,plusFace,tipF,options)
% Forward-ray intersections are in the DEFORMED frame. tipF maps local
% reference increments to deformed increments of the mean crack-face trace.
% Return reference-grid growth, capped by maxExtension at each tip.
if nargin<5
    options = struct;
end
o = struct('maxExtension',10,'tangentLength',10,'parallelTolerance',1e-6);
for name = fieldnames(options)'
    assert(isfield(o,name{1}),'Unknown face-growth option: %s',name{1});
    o.(name{1}) = options.(name{1});
end
for name = fieldnames(o)'
    validateattributes(o.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(path,{'numeric'},{'2d','ncols',2,'finite','real'});
assert(size(path,1)>=2 && isequal(size(path),size(minusFace),size(plusFace)), ...
    'Both traces must correspond to the ordered reference crack.');
assert(isequal(size(tipF),[2 2 2]),'Two tip deformation gradients are required.');
template = struct('accepted',false,'reason','no-intersection', ...
    'intersectionDeformed',[NaN NaN],'rayParameters',[NaN NaN], ...
    'targetReference',[NaN NaN],'extension',0,'uncappedExtension',0);
info = repmat(template,2,1);
targets = path([1 end],:);
for tip = 1:2
    ids = 1:size(path,1);
    if tip==2
        ids = fliplr(ids);
    end
    reference = path(ids,:);
    minus = minusFace(ids,:);
    plus = plusFace(ids,:);
    arc = [0;cumsum(vecnorm(diff(reference),2,2))];
    [arc,uniqueIDs] = unique(arc,'stable');
    if numel(arc)<2
        info(tip).reason = 'degenerate-path';
        continue
    end
    span = min(o.tangentLength,arc(end));
    inwardMinus = interp1(arc,minus(uniqueIDs,:),span);
    inwardPlus = interp1(arc,plus(uniqueIDs,:),span);
    dm = minus(1,:)-inwardMinus;
    dp = plus(1,:)-inwardPlus;
    F = tipF(:,:,tip);
    if any(~isfinite([dm dp F(:)'])) || min(norm(dm),norm(dp))<1e-9
        info(tip).reason = 'invalid-tangent';
        continue
    end
    dm = dm/norm(dm);
    dp = dp/norm(dp);
    A = [dm' -dp'];
    if abs(det(A))<=o.parallelTolerance
        info(tip).reason = 'parallel-rays';
        continue
    end
    parameters = A\(plus(1,:)-minus(1,:))';
    info(tip).rayParameters = parameters';
    if any(parameters<=1e-9)
        info(tip).reason = 'intersection-not-forward-on-both-faces';
        continue
    end
    intersection = minus(1,:)+parameters(1)*dm;
    info(tip).intersectionDeformed = intersection;
    if rcond(F)<1e-8 || det(F)<=0
        info(tip).reason = 'invalid-local-deformation';
        continue
    end
    midpoint = (minus(1,:)+plus(1,:))/2;
    increment = (F\(intersection-midpoint)')';
    inwardReference = interp1(arc,reference(uniqueIDs,:),span);
    outward = reference(1,:)-inwardReference;
    lengthHere = norm(increment);
    if lengthHere<=1e-9 || dot(increment,outward)<=0
        info(tip).reason = 'no-outward-reference-growth';
        continue
    end
    step = min(lengthHere,o.maxExtension);
    targets(tip,:) = reference(1,:)+step*increment/lengthHere;
    info(tip).accepted = true;
    info(tip).reason = 'forward-intersection';
    info(tip).targetReference = targets(tip,:);
    info(tip).extension = step;
    info(tip).uncappedExtension = lengthHere;
end
first = lineSegment(targets(1,:),path(1,:));
last = lineSegment(path(end,:),targets(2,:));
grown = [first(1:end-1,:);path;last(2:end,:)];
end

function points = lineSegment(a,b)
distance = norm(b-a);
if distance<1e-12
    points = a;
else
    points = a+linspace(0,1,ceil(distance)+1)'.*(b-a);
end
end
