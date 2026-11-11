function [warped,gx,gy,valid] = paperPFWarp(f,x,y,u)
%PAPERPFWARP Unweighted f(x-u), with exact derivatives of bilinear interpolation.
qx = x-u(:,1);
qy = y-u(:,2);
valid = qx>=1 & qx<=size(f,2) & qy>=1 & qy<=size(f,1);
warped = nan(size(x));
gx = warped;
gy = warped;
ix = min(floor(qx(valid)),size(f,2)-1);
iy = min(floor(qy(valid)),size(f,1)-1);
sx = qx(valid)-ix;
sy = qy(valid)-iy;
p00 = f(sub2ind(size(f),iy,ix));
p10 = f(sub2ind(size(f),iy,ix+1));
p01 = f(sub2ind(size(f),iy+1,ix));
p11 = f(sub2ind(size(f),iy+1,ix+1));
warped(valid) = (1-sx).*(1-sy).*p00+sx.*(1-sy).*p10+ ...
    (1-sx).*sy.*p01+sx.*sy.*p11;
gx(valid) = (1-sy).*(p10-p00)+sy.*(p11-p01);
gy(valid) = (1-sx).*(p01-p00)+sx.*(p11-p10);
end
