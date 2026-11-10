function [distance,segment,t,offset] = crackProjection2D(pathXY,pointsXY)
%CRACKPROJECTION2D Exact nearest-point projection onto a finite open polyline.
validateattributes(pathXY,{'numeric'},{'real','finite','2d','ncols',2});
validateattributes(pointsXY,{'numeric'},{'real','finite','2d','ncols',2});
pathXY=double(pathXY); pointsXY=double(pointsXY);
if size(pathXY,1)<2 || any(sum(diff(pathXY).^2,2)<1e-12)
    error('crackProjection2D:Path','Use at least two distinct consecutive points.');
end
n=size(pointsXY,1);
best=inf(n,1); segment=ones(n,1); t=zeros(n,1); offset=zeros(n,2);
for k=1:size(pathXY,1)-1
    tangent=pathXY(k+1,:)-pathXY(k,:);
    v=pointsXY-pathXY(k,:);
    tk=max(0,min(1,(v*tangent')/sum(tangent.^2)));
    delta=v-tk.*tangent;
    squared=sum(delta.^2,2);
    use=squared<best;
    best(use)=squared(use); segment(use)=k; t(use)=tk(use); offset(use,:)=delta(use,:);
end
distance=sqrt(best);
end
