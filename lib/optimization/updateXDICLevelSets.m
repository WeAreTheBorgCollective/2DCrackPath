function [updated,info] = updateXDICLevelSets(LS,R,mask,options)
%UPDATEXDICLEVELSETS Residual-driven normal propagation on the reference grid.
% R=f(X)-g(X+uXDIC); original images, no damage. Three signed-distance
% fields: phi (surface), psiStart and psiEnd (orthogonal endpoint fronts).
% Godunov propagation follows Rethore Eqs.25-32. Requested gradient velocity
% F=clip(gain*grad(mean(abs(R))).normal) is an adaptation of Eqs.35-36,
% whose literal normal-line maximum search is available as 'maximum'.
% Front positions follow the moving surface normally; no independent
% tangential growth law is assumed. All three fields are reconstructed from
% the accepted zero contour instead of using reinitialization PDE Eq.33.
% Full-window residual screening is optional; [1 1] uses raw |R|.
% This proposes geometry only; evolveXDIC re-solves and checks global RMS.
if nargin < 4
    options = struct();
end
o=struct('searchRadius',8,'sampleStep',0.25,'tangentSpan',5, ...
    'screeningWindowSize',[5 5],'velocityMethod','gradient', ...
    'gradientGain',0.1,'maxNormalStep',1,'smoothLength',3,'duration',1,'cfl',0.45,'maxBacktracks',4);
names=fieldnames(options);
for k=1:numel(names)
    if ~isfield(o,names{k})
        error('rethore:Option','Unknown option: %s',names{k});
    end
    o.(names{k})=options.(names{k});
end
for key={'searchRadius','sampleStep','tangentSpan','duration'}
    validateattributes(o.(key{1}),{'numeric'},{'scalar','finite','real','positive'});
end
validateattributes(o.smoothLength,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(o.cfl,{'numeric'},{'scalar','finite','real','positive','<=',0.5});
validateattributes(o.maxBacktracks,{'numeric'},{'scalar','integer','nonnegative','finite'});
if ~isequal(size(R),LS.imageSize,size(mask)) || ~islogical(mask) || ...
        any(~isfinite(R(mask))) || ~crackPathInDomain(LS.pathXY,mask)
    error('rethore:Domain','Finite raw residual and a valid initial crack inside a logical ROI are required.');
end
validateattributes(o.gradientGain,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(o.maxNormalStep,{'numeric'},{'scalar','real','finite','positive'});
o.velocityMethod=validatestring(o.velocityMethod,{'gradient','maximum'});
path=LS.pathXY;
s=[0;cumsum(vecnorm(diff(path),2,2))];
screenInput=double(R);
screenInput(~mask)=NaN;
[screened,~,~]=screenResidual2D(screenInput,o.screeningWindowSize,0,0,0);
screenMask=isfinite(screened);
% Invalid full-window centers contribute zero to the whole-path score.
% They never supply search targets or propagation speed; no raw-R fallback.
scoreMap=screened;
scoreMap(~screenMask)=0;
before=pathScore(path,scoreMap,mask);

raw=screened;
% NaNs deliberately invalidate finite-difference stencils at ROI boundaries.
[gradientX,gradientY]=gradient(screened);
gradientSupport=isfinite(gradientX) & isfinite(gradientY);
n=size(path,1);
supported=false(n,1);
offset=zeros(n,1);
targets=path;
normals=zeros(n,2);
sampleCount=max(1,ceil(o.searchRadius/o.sampleStep));
samples=linspace(-o.searchRadius,o.searchRadius,2*sampleCount+1)';
center=sampleCount+1;
for k=1:n
    lo=interp1(s,path,max(0,s(k)-o.tangentSpan));
    hi=interp1(s,path,min(s(end),s(k)+o.tangentSpan));
    tangent=(hi-lo)/max(norm(hi-lo),eps);
    normal=[-tangent(2),tangent(1)];
    normals(k,:)=normal;
    if strcmp(o.velocityMethod,'gradient')
        gx=interp2(gradientX,path(k,1),path(k,2),'linear',NaN);
        gy=interp2(gradientY,path(k,1),path(k,2),'linear',NaN);
        if ~isfinite(gx) || ~isfinite(gy)
            continue
        end
        supported(k)=true;
        normalGradient=gx*normal(1)+gy*normal(2);
        offset(k)=max(-o.maxNormalStep,min(o.maxNormalStep, ...
            o.gradientGain*normalGradient));
        targets(k,:)=path(k,:)+offset(k)*normal;
        continue
    end
    q=path(k,:)+samples.*normal;
    v=interp2(raw,q(:,1),q(:,2),'linear',NaN);
    allowed=interp2(double(mask),q(:,1),q(:,2),'linear',0)>1-1e-10 & isfinite(v);
    % Search only the contiguous valid interval containing the original point.
    if ~allowed(center)
        continue
    end
    supported(k)=true;
    a=center;
    b=center;
    while a>1 && allowed(a-1)
        a=a-1;
    end
    while b<numel(samples) && allowed(b+1)
        b=b+1;
    end
    indices=(a:b)';
    peak=max(v(indices));
    ties=indices(v(indices)>=peak-1e-10*max(1,peak));
    [~,j]=min(abs(samples(ties)));
    chosen=ties(j);
    % No drift on a flat profile.
    if peak<=v(center)+1e-10*max(1,peak)
        continue
    end
    targets(k,:)=q(chosen,:);
    offset(k)=interp2(LS.phi,q(chosen,1),q(chosen,2),'linear');
end
velocity=offset;
if o.smoothLength>0
    for k=1:n
        weights=max(0,1-abs(s-s(k))/o.smoothLength);
        velocity(k)=sum(weights.*offset)/sum(weights);
    end
end
velocity(~supported)=0; % Smoothing must not activate unsupported centers.
updated=LS;
info=struct('accepted',false,'reason','stationary','meanResidualBefore',before, ...
    'meanResidualAfter',before, ...
    'targetOffsets',offset,'velocity',velocity,'targetPoints',targets, ...
    'duration',0,'maxPathMotion',0,'substeps',0, ...
    'screeningWindowSize',o.screeningWindowSize,'validScreeningCenters',nnz(screenMask), ...
    'meanScreeningBefore',before,'meanScreeningAfter',before, ...
    'velocityMethod',o.velocityMethod);
if max(abs(velocity))<1e-10
    return
end
% Extend speed off the curve by nearest-segment interpolation in a narrow band.
band=screenMask & LS.distance<=max(o.searchRadius,o.maxNormalStep*o.duration)+3;
if strcmp(o.velocityMethod,'gradient')
    band=band & gradientSupport;
end
[y,x]=find(band);
[~,seg,t]=crackProjection2D(path,[x y]);
F=zeros(size(R));
F(band)=(1-t).*velocity(seg)+t.*velocity(seg+1);
info.reason='no-admissible-improving-contour';
for attempt=0:o.maxBacktracks
    duration=o.duration/2^attempt;
    [phi,propagation]=advectLevelSetGodunov(LS.phi,F,duration,o.cfl);
    predictedTips=path([1 end],:)+duration*velocity([1 end]).*normals([1 end],:);
    candidate=finiteContour(phi,predictedTips,path,mask, ...
        max(o.searchRadius,o.maxNormalStep*o.duration));
    if isempty(candidate)
        continue
    end
    score=pathScore(candidate,scoreMap,mask);
    if score<=before+1e-8*max(1,before)
        continue
    end
    updated=crackLevelSets2D(candidate,size(R));
    info.meanScreeningAfter=score;
    info.accepted=true;
    info.reason='accepted';
    info.meanResidualAfter=score;
    info.duration=duration;
    info.substeps=propagation.substeps;
    info.maxPathMotion=max(crackProjection2D(path,candidate));
    return
end
end

function score=pathScore(path,R,mask)
s=[0;cumsum(vecnorm(diff(path),2,2))];
q=interp1(s,path,linspace(0,s(end),max(2,ceil(s(end)*2)+1))');
values=interp2(abs(double(R)),q(:,1),q(:,2),'linear',NaN);
valid=interp2(double(mask),q(:,1),q(:,2),'linear',0)>1-1e-10;
if any(~valid | ~isfinite(values))
    score=-Inf;
else
    score=mean(values);
end
end

function best=finiteContour(phi,tips,oldPath,mask,radius)
% Cut one connected zero-contour at the projections of both predicted tips.
C=contourc(double(phi),[0 0]);
cursor=1;
best=zeros(0,2);
bestCost=Inf;
oldLength=sum(vecnorm(diff(oldPath),2,2));
while cursor<size(C,2)
    count=C(2,cursor);
    contour=C(:,cursor+1:cursor+count)';
    cursor=cursor+count+1;
    contour=compact(contour);
    if size(contour,1)<2
        continue
    end
    [distance,segment,~,delta]=crackProjection2D(contour,tips);
    if any(distance>radius+2)
        continue
    end
    q=tips-delta;
    if segment(1)>segment(2)
        contour=flipud(contour);
        [~,segment,~,delta]=crackProjection2D(contour,tips);
        q=tips-delta;
    end
    if segment(1)==segment(2)
        continue
    end
    candidate=compact([q(1,:);contour(segment(1)+1:segment(2),:);q(2,:)]);
    if size(candidate,1)<2
        continue
    end
    s=[0;cumsum(vecnorm(diff(candidate),2,2))];
    % Avoid a short unrelated contour or an excursion far away from the crack.
    if s(end)<0.8*oldLength || s(end)>1.3*oldLength
        continue
    end
    candidate=interp1(s,candidate,linspace(0,s(end),max(2,ceil(s(end))+1))');
    if ~crackPathInDomain(candidate,mask)
        continue
    end
    deviations=crackProjection2D(oldPath,candidate);
    if max(deviations)>radius+2
        continue
    end
    cost=mean(deviations)+sum(distance);
    if cost<bestCost
        bestCost=cost;
        best=candidate;
    end
end
end
function p=compact(p)
if isempty(p)
    return
end
keep=1;
for k=2:size(p,1)
    if norm(p(k,:)-p(keep(end),:))>=0.12
        keep(end+1)=k; %#ok<AGROW>
    end
end
p=p(keep,:);
end

