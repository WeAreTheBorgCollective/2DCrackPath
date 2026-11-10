function LS = crackLevelSets2D(pathXY,imageSize)
%CRACKLEVELSETS2D Crack signed distance plus two endpoint-front distances.
% Coordinates: x=column, y=row, unit pixel spacing. phi is signed distance
% to the polyline extended by endpoint tangent rays; positive on its left
% (in image coordinates). psiStart/psiEnd are signed distances to planes
% through each tip, with outward normals -tStart/+tEnd. The physical crack
% is phi=0 between the two tips. Gradients are orthogonal at the tips.
% Tip planes are LOCAL coordinates; for a strongly folded path they need
% not globally delimit the curve. distance is the exact unsigned distance
% to the FINITE polyline and is used for finite-tip damage construction.
% This extends Rethore et al. (2008), Eqs.20-22, to two crack tips.
validateattributes(imageSize,{'numeric'},{'vector','numel',2,'integer','positive'});
[x,y]=meshgrid(1:imageSize(2),1:imageSize(1)); points=[x(:),y(:)];
[distance,segment,~,delta]=crackProjection2D(pathXY,points);
tangents=diff(double(pathXY)); lengths=sqrt(sum(tangents.^2,2));
tangents=tangents./lengths;
tLocal=tangents(segment,:);
side=tLocal(:,1).*delta(:,2)-tLocal(:,2).*delta(:,1);
signSide=ones(size(side)); signSide(side<0)=-1;
phi=distance.*signSide; best=distance;
% Extend the first/last segment outwards for a well-defined signed surface.
for tip=[1 2]
    if tip==1
        origin=pathXY(1,:); tangent=tangents(1,:); direction=-1;
    else
        origin=pathXY(end,:); tangent=tangents(end,:); direction=1;
    end
    v=points-origin; along=v*tangent';
    outward=direction*along>=0;
    normalDistance=tangent(1)*v(:,2)-tangent(2)*v(:,1);
    use=outward & abs(normalDistance)<=best;
    phi(use)=normalDistance(use); best(use)=abs(normalDistance(use));
end
psiStart=-(points-pathXY(1,:))*tangents(1,:)';
psiEnd=(points-pathXY(end,:))*tangents(end,:)';
LS=struct('phi',reshape(phi,imageSize),'psiStart',reshape(psiStart,imageSize), ...
    'psiEnd',reshape(psiEnd,imageSize),'distance',reshape(distance,imageSize), ...
    'pathXY',double(pathXY),'imageSize',imageSize, ...
    'tipTangents',[tangents(1,:);tangents(end,:)]);
end
