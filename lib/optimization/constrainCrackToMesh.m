function [path,info] = constrainCrackToMesh(path,Mesh)
% Correct subpixel ROI/continuous-mesh boundary mismatch (at most one pixel).
tri = triangulation(double(Mesh.con(:,2:4)),double(Mesh.pos(:,1:2)));
info = struct('valid',true,'projectedPoints',0,'maxProjection',0);
outside = find(isnan(pointLocation(tri,path)));
if ~isempty(outside)
    edges = freeBoundary(tri);
    a = tri.Points(edges(:,1),:);
    b = tri.Points(edges(:,2),:);
    vector = b-a;
    length2 = sum(vector.^2,2);
    for k = outside'
        parameter = max(0,min(1,sum((path(k,:)-a).*vector,2)./length2));
        closest = a+parameter.*vector;
        [distance,j] = min(vecnorm(closest-path(k,:),2,2));
        info.maxProjection = max(info.maxProjection,distance);
        if distance>1
            info.valid = false;
            return
        end
        point = closest(j,:);
        % Tiny inward offset avoids floating-point ambiguity exactly on an edge.
        attached = edgeAttachments(tri,edges(j,:));
        element = attached{1}(1);
        center = mean(tri.Points(tri.ConnectivityList(element,:),:),1);
        inward = center-point;
        point = point+1e-8*inward/max(norm(inward),eps);
        path(k,:) = point;
        info.projectedPoints = info.projectedPoints+1;
    end
end
keep = [true;vecnorm(diff(path),2,2)>0.1];
path = path(keep,:);
if size(path,1)<2
    info.valid = false;
    return
end
% Check the segments as well as vertices, including notch boundaries.
for k = 1:size(path,1)-1
    count = max(1,ceil(4*norm(path(k+1,:)-path(k,:))));
    samples = path(k,:)+linspace(0,1,count+1)'.*(path(k+1,:)-path(k,:));
    if any(isnan(pointLocation(tri,samples)))
        info.valid = false;
        return
    end
end
end
