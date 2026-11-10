function ok = crackPathInDomain(pathXY,mask)
%CRACKPATHINDOMAIN Reject collapsed, crossing or out-of-domain trial curves.
ok=false;
if size(pathXY,1)<2 || size(pathXY,2)~=2 || any(~isfinite(pathXY(:))), return; end
edges=diff(pathXY);
if any(sum(edges.^2,2)<0.01), return; end
for k=1:size(edges,1)
    n=max(1,ceil(2*norm(edges(k,:))));
    points=pathXY(k,:)+(0:n)'/n.*edges(k,:);
    xy=round(points);
    if any(xy(:,1)<1 | xy(:,1)>size(mask,2) | xy(:,2)<1 | xy(:,2)>size(mask,1)), return; end
    if ~all(mask(sub2ind(size(mask),xy(:,2),xy(:,1)))), return; end
end
% Non-neighboring segment intersections, including collinear overlap.
cross2=@(a,b) a(1)*b(2)-a(2)*b(1);
for i=1:size(edges,1)-2
    for j=i+2:size(edges,1)
        a=pathXY(i,:); b=pathXY(j,:); u=edges(i,:); v=edges(j,:);
        den=cross2(u,v); delta=b-a;
        if abs(den)>1e-12
            s=cross2(delta,v)/den; t=cross2(delta,u)/den;
            if s>=0 && s<=1 && t>=0 && t<=1, return; end
        elseif abs(cross2(delta,u))<1e-10
            span=sort([dot(delta,u),dot(delta+v,u)])/dot(u,u);
            if span(1)<=1 && span(2)>=0, return; end
        end
    end
end
ok=true;
end
