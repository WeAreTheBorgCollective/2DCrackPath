function model = paperPFModel(mask,spacing)
%PAPERPFMODEL Coarse Q1 displacement and pixel-spaced lumped T3 damage FE.
% Both meshes live in the DEFORMED x frame used by the paper.
validateattributes(spacing,{'numeric'},{'scalar','integer','>=',1});
[y,x] = find(mask);
if isempty(x)
    error('paperPF:ROI','The ROI is empty.');
end
model.mask = mask;
model.ids = find(mask);
model.x = double(x);
model.y = double(y);
model.colors = mod(x+y,2); % Two independent sets of the unit-grid stiffness.
model.imageSize = size(mask);
xn = min(x):spacing:max(x);
yn = min(y):spacing:max(y);
if xn(end)~=max(x)
    xn(end+1) = max(x);
end
if yn(end)~=max(y)
    yn(end+1) = max(y);
end
if numel(xn)<2 || numel(yn)<2
    error('paperPF:ROI','ROI must contain at least two rows and columns.');
end
[X,Y] = meshgrid(xn,yn);
model.nodes = [X(:) Y(:)];
ix = min(floor((x-xn(1))/spacing)+1,numel(xn)-1);
iy = min(floor((y-yn(1))/spacing)+1,numel(yn)-1);
hx = reshape(xn(ix+1)-xn(ix),[],1);
hy = reshape(yn(iy+1)-yn(iy),[],1);
sx = (x-reshape(xn(ix),[],1))./hx;
sy = (y-reshape(yn(iy),[],1))./hy;
n1 = sub2ind(size(X),iy,ix);
nodes = [n1 n1+size(X,1) n1+size(X,1)+1 n1+1];
N = [(1-sx).*(1-sy) sx.*(1-sy) sx.*sy (1-sx).*sy];
Bx = [-(1-sy)./hx (1-sy)./hx sy./hx -sy./hx];
By = [-(1-sx)./hy -sx./hy sx./hy (1-sx)./hy];
rows = repmat((1:numel(x))',1,4);
model.N = sparse(rows(:),nodes(:),N(:),numel(x),numel(X));
model.Bx = sparse(rows(:),nodes(:),Bx(:),numel(x),numel(X));
model.By = sparse(rows(:),nodes(:),By(:),numel(x),numel(X));
active = full(sum(model.N,1))>0;
model.N = model.N(:,active);
model.Bx = model.Bx(:,active);
model.By = model.By(:,active);
model.nodes = model.nodes(active,:);
% Two unit right triangles per complete pixel cell, diagonal lower-left to
% upper-right in matrix indexing. Exact T3 stiffness; vertex-lumped mass.
id = zeros(size(mask));
id(mask) = 1:nnz(mask);
cells = mask(1:end-1,1:end-1) & mask(2:end,1:end-1) & ...
    mask(1:end-1,2:end) & mask(2:end,2:end);
[cy,cx] = find(cells);
v1 = id(sub2ind(size(mask),cy,cx));
v2 = id(sub2ind(size(mask),cy,cx+1));
v3 = id(sub2ind(size(mask),cy+1,cx+1));
v4 = id(sub2ind(size(mask),cy+1,cx));
triangles = [v1 v2 v3;v1 v3 v4];
model.mass = accumarray(triangles(:),1/6,[nnz(mask) 1]);
if any(model.mass==0)
    error('paperPF:ROI','Every ROI pixel must belong to a complete unit cell.');
end
% Each cell contributes half-stiffness on each of its four outer edges.
edgeA = [v1;v2;v3;v4];
edgeB = [v2;v3;v4;v1];
weights = 0.5*ones(size(edgeA));
model.K = sparse([edgeA;edgeB;edgeA;edgeB], ...
    [edgeA;edgeB;edgeB;edgeA],[weights;weights;-weights;-weights],nnz(mask),nnz(mask));
end
