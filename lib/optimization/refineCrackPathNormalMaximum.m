function [refined,G,pathMask,info] = refineCrackPathNormalMaximum(path,R,options,Mesh)
% Relocate every ordered crack-path node to the maximum |R| on its normal.
% All normals use the ORIGINAL path; updates are simultaneous.
if nargin<3
    options = struct;
end
o = struct('radius',30,'density',1);
for name = fieldnames(options)'
    assert(isfield(o,name{1}),'Unknown normal-search option: %s',name{1});
    o.(name{1}) = options.(name{1});
end
validateattributes(o.radius,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(o.density,{'numeric'},{'scalar','finite','positive'});
validateattributes(path,{'numeric'},{'2d','ncols',2,'finite','real','nonempty'});
assert(size(path,1)>=2,'At least two crack nodes are needed.');
validateattributes(R,{'numeric'},{'2d','real','nonempty'});
count = floor(o.radius*o.density+1e-12);
offsets = (-count:count)/o.density;
mask = isfinite(R);
values = abs(double(R));
n = size(path,1);
refined = double(path);
normals = zeros(n,2);
scores = nan(n,numel(offsets));
selectedOffsets = zeros(n,1);
before = sampleROIImage(values,mask,path);
after = before;
if nargin>=4 && ~isempty(Mesh)
    tri = triangulation(double(Mesh.con(:,2:4)),double(Mesh.pos(:,1:2)));
else
    tri = [];
end
for k = 1:n
    % Centered chord between adjacent distinct nodes; one-sided at endpoints.
    left = find(vecnorm(path(1:k-1,:)-path(k,:),2,2)>1e-12,1,'last');
    right = find(vecnorm(path(k+1:end,:)-path(k,:),2,2)>1e-12,1);
    if isempty(left)
        left = k;
    end
    if isempty(right)
        right = k;
    else
        right = k+right;
    end
    tangent = path(right,:)-path(left,:);
    if norm(tangent)<1e-12
        continue
    end
    tangent = tangent/norm(tangent);
    normal = [-tangent(2) tangent(1)];
    normals(k,:) = normal;
    points = path(k,:)+offsets'.*normal;
    v = sampleROIImage(values,mask,points);
    if ~isempty(tri)
        v(isnan(pointLocation(tri,points))) = NaN;
    end
    scores(k,:) = v';
    valid = find(isfinite(v));
    if isempty(valid)
        continue
    end
    peak = max(v(valid));
    ties = valid(v(valid)==peak);
    [~,j] = min(abs(offsets(ties)));
    selected = ties(j); % nearest maximum; zero wins on a flat profile
    refined(k,:) = points(selected,:);
    selectedOffsets(k) = offsets(selected);
    after(k) = v(selected);
end
weights = vecnorm(diff(refined),2,2);
G = graph((1:n-1)',(2:n)',weights,n);
G.Nodes.X = refined(:,1);
G.Nodes.Y = refined(:,2);
G.Nodes.OriginalX = path(:,1);
G.Nodes.OriginalY = path(:,2);
G.Nodes.NormalOffset = selectedOffsets;
pathMask = false(size(R));
for k = 1:n-1
    steps = max(1,ceil(max(abs(refined(k+1,:)-refined(k,:)))));
    pixels = round(refined(k,:)+linspace(0,1,steps+1)'.*(refined(k+1,:)-refined(k,:)));
    inside = pixels(:,1)>=1 & pixels(:,1)<=size(R,2) & ...
        pixels(:,2)>=1 & pixels(:,2)<=size(R,1);
    pixels = pixels(inside,:);
    pathMask(sub2ind(size(R),pixels(:,2),pixels(:,1))) = true;
end
info = struct('options',o,'offsets',offsets,'normals',normals, ...
    'sampleValues',scores,'selectedOffsets',selectedOffsets, ...
    'residualBefore',before,'residualAfter',after, ...
    'pathLength',sum(weights),'validForLevelSets',crackPathInDomain(refined,mask), ...
    'scoreSource','absolute raw residual, ROI-weighted bilinear samples');
end
