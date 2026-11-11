function prep = prepareXDIC(Mesh,imageSize,crackPathXY)
%PREPAREXDIC Pixel basis for 2D T3 FE with shifted-Heaviside enrichment.
% crackPathXY is in REFERENCE image coordinates [column,row].
% Tip-containing nodal supports are excluded; there are no tip functions.
if nargin<3 || isempty(crackPathXY)
    crackPathXY = zeros(0,2);
end
if exist('xdicPrepMex','file')~=3
    error('xdic:Build','Run buildXDICKernels first.');
end
validateattributes(Mesh.pos,{'numeric'},{'real','finite','2d'});
validateattributes(Mesh.con,{'numeric'},{'real','finite','2d'});
if size(Mesh.pos,2)<2 || size(Mesh.pos,2)>3 || ...
        (size(Mesh.pos,2)==3 && any(abs(Mesh.pos(:,3)-Mesh.pos(1,3))>1e-10))
    error('xdic:Mesh','Only planar XY meshes are supported.');
end
if size(Mesh.con,2)~=4 || any(Mesh.con(:,1)~=0)
    error('xdic:Mesh','Use Correli T3 connectivity [0 node1 node2 node3].');
end
validateattributes(imageSize,{'numeric'},{'vector','numel',2,'integer','positive','finite'});
validateattributes(crackPathXY,{'numeric'},{'real','finite','2d','ncols',2});
prep = xdicPrepMex(double(Mesh.pos(:,1:2)),double(Mesh.con(:,2:4)), ...
    double(imageSize(:)'),double(crackPathXY));
end
