function Mesh = loadScreeningMesh(filename,imageSize,offset)
% Convert MAT mesh.pos/con to Correli T3 layout, preserving image coordinates.
[~,~,extension] = fileparts(filename);
if strcmpi(extension,'.mat')
    data = load(filename);
    if isfield(data,'mesh')
        source = data.mesh;
    elseif isfield(data,'Mesh')
        source = data.Mesh;
    else
        error('screening:Mesh','Expected mesh or Mesh in %s.',filename);
    end
else
    Mesh = struct;
    run(filename);
    source = Mesh;
end
assert(size(source.pos,2)>=2 && size(source.con,2)==4, ...
    'Expected XY[Z] positions and [elementType n1 n2 n3] connectivity.');
xy = double(source.pos(:,1:2));
tri = double(source.con(:,2:4));
assert(all(isfinite(xy(:))) && all(isfinite(tri(:))) && ...
    all(tri(:)==fix(tri(:))) && all(tri(:)>=1) && ...
    all(tri(:)<=size(xy,1)),'Invalid mesh positions or node indices.');
if ischar(offset) || isstring(offset)
    assert(strcmp(offset,'original'),'Unknown mesh offset.');
    offset = [imageSize(2)/2 imageSize(1)/2-90];
end
xy = xy + offset;
v1 = xy(tri(:,2),:)-xy(tri(:,1),:);
v2 = xy(tri(:,3),:)-xy(tri(:,1),:);
area2 = v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);
assert(all(abs(area2)>eps),'Degenerate triangle found.');
% Correli requires counterclockwise triangles.
reverse = area2<0;
tri(reverse,[2 3]) = tri(reverse,[3 2]);
Mesh.pos = [xy zeros(size(xy,1),1)];
Mesh.con = [zeros(size(tri,1),1) tri];
end
