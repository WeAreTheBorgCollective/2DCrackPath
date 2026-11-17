function [extended,mask,info] = extendCrackPathToMask(path,binary,tangentLength)
% Extend each ordered-path tip to the furthest forward foreground pixel.
% Candidates belong to the tip's 8-connected binary component and lie within
% one pixel of its outward ray. Connections are straight, including over holes.
if nargin<3
    tangentLength = 10;
end
validateattributes(path,{'numeric'},{'2d','ncols',2,'finite','real'});
validateattributes(binary,{'logical'},{'2d'});
validateattributes(tangentLength,{'numeric'},{'scalar','positive','finite'});
extended = double(path);
mask = false(size(binary));
info = struct('tangentLength',tangentLength,'lineTolerance',1, ...
    'originalEndpoints',zeros(0,2),'endpoints',zeros(0,2), ...
    'directions',zeros(2,2),'extensionLengths',[0;0], ...
    'componentIDs',[0;0],'reason',{{'insufficient-path','insufficient-path'}});
if size(path,1)<2
    return
end
assert(all(path(:,1)>=1 & path(:,1)<=size(binary,2) & ...
    path(:,2)>=1 & path(:,2)<=size(binary,1)),'Path is outside the image.');
labels = bwlabel(binary,8);
info.originalEndpoints = path([1 end],:);
info.endpoints = info.originalEndpoints;
for tip = 1:2
    ordered = double(path);
    if tip==2
        ordered = flipud(ordered);
    end
    endpoint = ordered(1,:);
    distance = [0;cumsum(vecnorm(diff(ordered),2,2))];
    keep = [true;diff(distance)>0];
    ordered = ordered(keep,:);
    distance = distance(keep);
    if numel(distance)<2
        info.reason{tip} = 'zero-length-path';
        continue
    end
    interior = interp1(distance,ordered,min(tangentLength,distance(end)));
    direction = endpoint-interior;
    if norm(direction)<=eps
        info.reason{tip} = 'undefined-direction';
        continue
    end
    direction = direction/norm(direction);
    info.directions(tip,:) = direction;
    pixel = round(endpoint);
    component = labels(pixel(2),pixel(1));
    info.componentIDs(tip) = component;
    if component==0
        info.reason{tip} = 'tip-outside-binary-mask';
        continue
    end
    [y,x] = find(labels==component);
    delta = [x y]-endpoint;
    forward = delta*direction';
    perpendicular = abs(delta(:,1)*direction(2)-delta(:,2)*direction(1));
    eligible = find(forward>1e-9 & perpendicular<=1+1e-9);
    if isempty(eligible)
        info.reason{tip} = 'no-forward-pixel';
        continue
    end
    [~,order] = sortrows([-forward(eligible) perpendicular(eligible)],[1 2]);
    selected = eligible(order(1));
    target = [x(selected) y(selected)];
    info.endpoints(tip,:) = target;
    info.extensionLengths(tip) = norm(target-endpoint);
    info.reason{tip} = 'extended';
end
% Evaluate both directions on the original path, then append both extensions.
startSegment = segment(info.endpoints(1,:),path(1,:));
endSegment = segment(path(end,:),info.endpoints(2,:));
extended = [startSegment(1:end-1,:);path;endSegment(2:end,:)];
for k = 1:size(extended,1)-1
    pixels = round(segment(extended(k,:),extended(k+1,:)));
    mask(sub2ind(size(mask),pixels(:,2),pixels(:,1))) = true;
end
info.originalLength = sum(vecnorm(diff(path),2,2));
info.extendedLength = sum(vecnorm(diff(extended),2,2));
end

function points = segment(a,b)
if norm(b-a)<=eps
    points = a;
    return
end
steps = max(1,ceil(norm(b-a)));
points = a+linspace(0,1,steps+1)'.*(b-a);
end
