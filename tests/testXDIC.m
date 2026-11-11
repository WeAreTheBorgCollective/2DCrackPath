function testXDIC
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setupPaths;
[x,y] = meshgrid(8:8:72);
Mesh.pos = [x(:),y(:),zeros(numel(x),1)];
tri = delaunay(x(:),y(:));
Mesh.con = [zeros(size(tri,1),1),tri];
p = prepareXDIC(Mesh,[80 80],[]);
assert(isempty(p.enrichedNodes));
assert(max(abs(sum(p.basis,2)-1))<1e-12);
[py,px] = ind2sub([80 80],p.pixelIndex);
assert(max(abs(p.basis*Mesh.pos(:,1)-px))<1e-12);
crack = [35.5 1;35.5 80];
q = prepareXDIC(Mesh,[80 80],crack);
assert(~isempty(q.enrichedNodes));
qr = prepareXDIC(Mesh,[80 80],flipud(crack));
nn = size(Mesh.pos,1);
assert(norm(q.basis(:,nn+1:end)+qr.basis(:,nn+1:end),'fro')<1e-12);
qt = prepareXDIC(Mesh,[80 80],[35.5 20;35.5 60]);
shared = intersect(qt.enrichedNodes,qt.tipSupportNodes);
assert(~isempty(shared));
fullNodes = unique(Mesh.con(qt.fullyCutElements,2:4));
assert(all(ismember(qt.enrichedNodes,fullNodes)));
tiny = prepareXDIC(Mesh,[80 80],[34 34;35 34.5]);
assert(isempty(tiny.enrichedNodes));
[X,Y] = meshgrid(1:80);
texture = @(x,y) 100+25*sin(.23*x+.13*y)+20*cos(.17*x-.27*y)+15*sin(.37*y);
f = texture(X,Y);
g = texture(X,Y-.6*(X>35.5));
s = xdicCorrMex(q,f,g,zeros(2*size(q.basis,2),1),'linear',true);
gx = (f(q.pixelIndex+80)-f(q.pixelIndex-80))/2;
gy = (f(q.pixelIndex+1)-f(q.pixelIndex-1))/2;
J = zeros(numel(q.pixelIndex),2*size(q.basis,2));
J(:,1:2:end) = q.basis.*gx;
J(:,2:2:end) = q.basis.*gy;
assert(norm(s.H-J'*J,'fro')/norm(s.H,'fro')<1e-12);
assert(norm(s.b-J'*s.R(q.pixelIndex))/norm(s.b)<1e-12);
opts = struct('maxIterations',60);
[~,rFE] = correlateXDIC(f,g,Mesh,[],[],opts);
[u,rXF] = correlateXDIC(f,g,Mesh,crack,[],opts);
uy = u.U(:,:,2);
mask = isfinite(uy) & abs(X-35.5)>5;
errorUY = sqrt(mean((uy(mask)-.6*(X(mask)>35.5)).^2));
assert(rXF<.2*rFE);
assert(errorUY<.03);
fprintf('PASS: basis, tip exclusion, orientation, full Hessian/RHS, jump recovery. FE %.6g XF %.6g displacement error %.6g px\n',rFE,rXF,errorUY);
% Compare against the installed Correli kernel on one noninteger triangle.
if exist('kernel_corr','file')
    M.pos = [10.2 11.3 0;65.6 13.1 0;20.4 67.2 0];
    M.con = [0 1 2 3];
    p = prepareXDIC(M,size(f),[]);
    s = xdicCorrMex(p,single(f),single(g),zeros(6,1),'linear',true);
    native = kernel_corr('linear',kernel_prep(M,size(f)));
    [H,b,R] = native(single(f),single(g),zeros(9,1),'MbR');
    ix = [1 2 4 5 7 8];
    common = isfinite(R) & isfinite(s.R);
    herror = norm(H(ix,ix)-s.H,'fro')/norm(s.H,'fro');
    berror = norm(b(ix)-s.b)/norm(s.b);
    rerror = max(abs(double(R(common))-s.R(common)));
    fprintf('Native parity: Hessian %.3g RHS %.3g residual %.3g\n',herror,berror,rerror);
    assert(herror<1e-5 && berror<1e-5 && rerror<1e-4);
end
end


