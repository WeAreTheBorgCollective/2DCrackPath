function [pathXY, pathMask, G, info] = connectCrackGraph(skeleton, options)
%CONNECTCRACKGRAPH Euclidean pixel graph, gap merging and one principal path.
% Strini & Schiavi (2022), Sections 3.2, 3.3.1, 3.3.2 and 3.3.5.
% Nodes are skeleton pixels; 8-neighbor edges have Euclidean weights.
% Small/compact components are filtered, then the nearest available pair
% of endpoints in different components is joined, iteratively, without
% cycles, up to MaxGap. Only segments entirely in ValidMask may be joined.
% This conservative merging policy recomputes endpoints after each join.
% The longest endpoint-to-endpoint weighted shortest path is returned.
% If several components remain, the longest such path among them wins;
% disconnected components/side branches are NOT silently claimed as joined.
% Components without two endpoints use all their vertices as path terminals.
% This is a principal-path adaptation, not the paper's full refinement
% pipeline: subgraph recovery, secondary branches and correlation omitted.
%
% options.MaxGap = 100 [px]; Inf explicitly allows arbitrarily long bridges.
% options.MinComponentDiameter = 25 [px] (endpoint Euclidean separation).
% options.MaxVertexDensity = 2 [vertices/px]; Inf disables density filtering.
% options.ValidMask = true(size(skeleton)).
% pathXY: ordered [x=column,y=row] nodes, including distal-edge endpoints.
% pathMask: rasterized polyline, including interpolated gaps.
% G: retained/merged graph, Nodes.X/Y/PixelIndex, Edges.Weight/IsBridge.
% info.pathLength counts true Euclidean edge lengths, not raster pixels.
if nargin < 2, options = struct(); end
validateattributes(skeleton,{'logical'},{'2d','nonempty'});
defaults = struct('MaxGap',100,'MinComponentDiameter',25, ...
    'MaxVertexDensity',2,'ValidMask',true(size(skeleton)));
names = fieldnames(options);
for k=1:numel(names)
    if ~isfield(defaults,names{k})
        error('connectCrackGraph:Option','Unknown option: %s',names{k});
    end
    defaults.(names{k}) = options.(names{k});
end
o = defaults;
validateattributes(o.MaxGap,{'numeric'},{'scalar','real','nonnegative','nonnan'});
validateattributes(o.MinComponentDiameter,{'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(o.MaxVertexDensity,{'numeric'},{'scalar','real','positive','nonnan'});
if ~islogical(o.ValidMask) || ~isequal(size(o.ValidMask),size(skeleton))
    error('connectCrackGraph:Mask','ValidMask must be logical and match the skeleton.');
end
if any(skeleton(:) & ~o.ValidMask(:))
    error('connectCrackGraph:Mask','Skeleton pixels must lie inside ValidMask.');
end
G = pixelGraph(skeleton);
info = struct('initialComponents',0,'removedComponents',0, ...
    'remainingComponents',0,'bridgeCount',0,'pathLength',0, ...
    'pathNodeIDs',zeros(0,1),'excludedNodeCount',0,'allRetainedConnected',false);
pathXY = zeros(0,2); pathMask = false(size(skeleton));
if numnodes(G)==0, return; end

labels = conncomp(G);
info.initialComponents = max(labels);
remove = false(numnodes(G),1);
for c=1:max(labels)
    ids = find(labels==c);
    ends = ids(degree(G,ids)<=1);
    if numel(ends)<2, ends=ids; end
    diameter = 0;
    for j=1:numel(ends)
        diameter = max(diameter,max(hypot( ...
            G.Nodes.X(ends)-G.Nodes.X(ends(j)), ...
            G.Nodes.Y(ends)-G.Nodes.Y(ends(j)))));
    end
    density = numel(ids)/max(diameter,eps);
    if diameter<o.MinComponentDiameter || density>o.MaxVertexDensity
        remove(ids)=true;
        info.removedComponents=info.removedComponents+1;
    end
end
G = rmnode(G,find(remove));
if numnodes(G)==0, return; end

% Nearest feasible endpoint pair. Distal links are flagged, not confused
% with observed skeleton pixels. Recompute components/degrees after joining.
while true
    labels = conncomp(G);
    if max(labels)==1, break; end
    ends = find(degree(G)<=1);
    best = Inf; pair = [];
    for i=1:numel(ends)
        a = ends(i);
        for j=i+1:numel(ends)
            b = ends(j);
            if labels(a)==labels(b), continue; end
            d = hypot(G.Nodes.X(a)-G.Nodes.X(b),G.Nodes.Y(a)-G.Nodes.Y(b));
            if d>o.MaxGap || d>=best, continue; end
            xy = rasterLine([G.Nodes.X(a),G.Nodes.Y(a)], ...
                [G.Nodes.X(b),G.Nodes.Y(b)]);
            lin = sub2ind(size(skeleton),xy(:,2),xy(:,1));
            if all(o.ValidMask(lin))
                best=d; pair=[a b];
            end
        end
    end
    if isempty(pair), break; end
    edge = table(pair,best,true,'VariableNames',{'EndNodes','Weight','IsBridge'});
    G = addedge(G,edge);
end

labels = conncomp(G);
info.remainingComponents = max(labels);
info.allRetainedConnected = info.remainingComponents==1;
info.bridgeCount = nnz(G.Edges.IsBridge);
bestLength = -Inf; bestPath = [];
% Exact terminal geodesic diameter, one distance row at a time to avoid an
% all-node distance matrix. Loops use all nodes as the documented fallback.
for c=1:max(labels)
    ids=find(labels==c); terminals=ids(degree(G,ids)<=1);
    if numel(terminals)<2, terminals=ids; end
    for j=1:numel(terminals)
        d=distances(G,terminals(j),terminals,'Method','positive');
        [lengthHere,k]=max(d);
        if lengthHere>bestLength
            bestLength=lengthHere;
            bestPath=shortestpath(G,terminals(j),terminals(k),'Method','positive');
        end
    end
end
pathXY = [G.Nodes.X(bestPath),G.Nodes.Y(bestPath)];
% Stable orientation: smaller x first, breaking vertical ties by y.
if size(pathXY,1)>1 && (pathXY(1,1)>pathXY(end,1) || ...
        (pathXY(1,1)==pathXY(end,1) && pathXY(1,2)>pathXY(end,2)))
    bestPath=fliplr(bestPath); pathXY=flipud(pathXY);
end
if size(pathXY,1)==1
    pathMask(pathXY(1,2),pathXY(1,1))=true;
end
for k=1:size(pathXY,1)-1
    xy=rasterLine(pathXY(k,:),pathXY(k+1,:));
    pathMask(sub2ind(size(skeleton),xy(:,2),xy(:,1)))=true;
end
info.pathNodeIDs=bestPath(:);
info.pathLength=max(0,bestLength);
info.excludedNodeCount=numnodes(G)-numel(bestPath);
end

function G=pixelGraph(BW)
[y,x]=find(BW); n=numel(x);
map=zeros(size(BW)); map(BW)=1:n;
s=zeros(0,1); t=s;
offsets=[0 1;1 -1;1 0;1 1];
for k=1:4
    yy=y+offsets(k,1); xx=x+offsets(k,2);
    valid=yy>=1 & yy<=size(BW,1) & xx>=1 & xx<=size(BW,2);
    a=find(valid); b=map(sub2ind(size(BW),yy(valid),xx(valid)));
    a=a(b>0); b=b(b>0);
    s=[s;a]; t=[t;b]; %#ok<AGROW>
end
G=graph(s,t,hypot(x(s)-x(t),y(s)-y(t)),n);
G.Nodes.X=x; G.Nodes.Y=y; G.Nodes.PixelIndex=find(BW);
G.Edges.IsBridge=false(numedges(G),1);
end

function xy=rasterLine(a,b)
steps=max(abs(b-a));
if steps==0, xy=a; return; end
fraction=(0:steps)'/steps;
xy=round(a+fraction.*(b-a));
end
