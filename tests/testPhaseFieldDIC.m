function tests=testPhaseFieldDIC
tests=functiontests(localfunctions);
end
function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'lib','optimization'),fullfile(root,'lib','corr'));
end
function testStraightLevelSetsAndTipOrthogonality(tc)
LS=crackLevelSets2D([10 20;40 20],[50 60]);
[x,y]=meshgrid(1:60,1:50);
verifyEqual(tc,LS.phi,y-20,'AbsTol',1e-12);
verifyEqual(tc,LS.psiStart,10-x,'AbsTol',1e-12);
verifyEqual(tc,LS.psiEnd,x-40,'AbsTol',1e-12);
verifyEqual(tc,LS.distance(20,5),5);
verifyEqual(tc,LS.distance(24,7),5);
[gx,gy]=gradient(LS.phi); [sx,sy]=gradient(LS.psiStart);
verifyEqual(tc,hypot(gx,gy),ones(size(x)),'AbsTol',1e-12);
verifyEqual(tc,gx.*sx+gy.*sy,zeros(size(x)),'AbsTol',1e-12);
end
function testCurvedSubpixelProjection(tc)
path=[3.5 5.5;8.5 5.5;8.5 10.5];
[r,seg,t]=crackProjection2D(path,[6 6.5;9.5 8;2.5 5.5]);
verifyEqual(tc,r,[1;1;1],'AbsTol',1e-12);
verifyEqual(tc,seg,[1;2;1]); verifyEqual(tc,t,[0.5;0.5;0],'AbsTol',1e-12);
LS=crackLevelSets2D(path,[15 15]);
verifyGreaterThan(tc,LS.phi(7,6),0);
verifyLessThan(tc,LS.phi(4,6),0);
end
function testDamageAndOutsidePreservation(tc)
LS=crackLevelSets2D([10 20;40 20],[50 60]);
mask=false(50,60); mask(5:45,3:55)=true;
base=nan(50,60); p=phaseFieldParameters;
d=createDamageField(LS,mask,p,base);
verifyTrue(tc,all(isnan(d(~mask)))); verifyEqual(tc,d(20,10:40),ones(1,31));
r=sqrt(5^2+p.distanceSmoothing^2)-p.distanceSmoothing;
verifyEqual(tc,d(25,25),exp(-sqrt(p.mu)*r/p.lc),'AbsTol',1e-12);
verifyGreaterThanOrEqual(tc,min(d(mask)),0); verifyLessThanOrEqual(tc,max(d(mask)),1);
file=[tempname '.mat']; clean=onCleanup(@()delete(file)); %#ok<NASGU>
createDamageField(LS,mask,p,base,file);
stored=load(file); verifyEqual(tc,stored.Damage,d);
verifyEqual(tc,stored.measurementMask,mask); verifyEqual(tc,stored.levelSets.pathXY,LS.pathXY);
end
function testEnergyDerivativeAndNoFlux(tc)
mask=true(7,9); mask(3,4)=false; mask(:,1)=false;
d=reshape(linspace(0.2,0.7,63),7,9);
R=reshape(cos(1:63),7,9); R(~mask)=NaN; d(~mask)=NaN;
p=phaseFieldParameters(struct('grayScale',3));
[E,gd,parts]=phaseFieldEnergy2D(d,R,mask,p);
verifyEqual(tc,E,parts.dic+parts.crack,'AbsTol',1e-12);
for id=[10 20 40 60]
    if ~mask(id), continue; end
    dp=d; dm=d; h=1e-6; dp(id)=dp(id)+h; dm(id)=dm(id)-h;
    numerical=(phaseFieldEnergy2D(dp,R,mask,p)-phaseFieldEnergy2D(dm,R,mask,p))/(2*h);
    verifyEqual(tc,gd(id),numerical,'AbsTol',1e-7);
end
constant=0.3*ones(size(mask));
[~,~,parts]=phaseFieldEnergy2D(constant,R,mask,p);
verifyEqual(tc,parts.gradientDamage,0);
verifyEqual(tc,gd(~mask),zeros(nnz(~mask),1));
end
function testShapeDerivativeAgainstFiniteDifference(tc)
path=[7.23 12.31;23.71 13.13]; dims=[30 32]; mask=true(dims);
[x,y]=meshgrid(1:dims(2),1:dims(1));
R=1+0.01*x+0.02*y; p=phaseFieldParameters(struct('maxBacktracks',0));
LS=crackLevelSets2D(path,dims);
[~,~,info]=updateCrackLevelSets(LS,R,mask,p);
h=1e-5;
for k=1:2
    for axis=1:2
        plus=path; minus=path; plus(k,axis)=plus(k,axis)+h; minus(k,axis)=minus(k,axis)-h;
        dp=createDamageField(crackLevelSets2D(plus,dims),mask,p);
        dm=createDamageField(crackLevelSets2D(minus,dims),mask,p);
        numerical=(phaseFieldEnergy2D(dp,R,mask,p)-phaseFieldEnergy2D(dm,R,mask,p))/(2*h);
        verifyEqual(tc,info.shapeGradient(k,axis),numerical,'AbsTol',1e-5);
    end
end
end
function testShapeUpdateAndTipRebuild(tc)
dims=[40 40]; mask=true(dims); [~,y]=meshgrid(1:40,1:40);
R=10*exp(-(y-22).^2/6);
LS=crackLevelSets2D([8 18;20 18;32 18],dims); p=phaseFieldParameters;
[new,d,info]=updateCrackLevelSets(LS,R,mask,p);
verifyTrue(tc,info.accepted); verifyLessThan(tc,info.energyAfter,info.energyBefore);
verifyGreaterThan(tc,new.pathXY(2,2),LS.pathXY(2,2));
verifyEqual(tc,d,createDamageField(new,mask,p),'AbsTol',1e-12);
verifyGreaterThan(tc,norm(new.psiStart-LS.psiStart,'fro'),0);
verifyGreaterThan(tc,norm(new.psiEnd-LS.psiEnd,'fro'),0);
verifyTrue(tc,crackPathInDomain(new.pathXY,mask));
verifyLessThanOrEqual(tc,max(vecnorm(new.pathXY-LS.pathXY,2,2)),p.maxGeometryStep+1e-12);
end
function testImageWeighting(tc)
f=single(reshape(1:30,5,6)); g=2*f;
mask=false(5,6); mask(2:4,2:5)=true;
d=0.4*ones(5,6); d(~mask)=NaN;
[fp,gp]=applyDamageToImages(f,g,d,mask);
verifyEqual(tc,fp(~mask),f(~mask)); verifyEqual(tc,gp(~mask),g(~mask));
verifyEqual(tc,(fp(mask)-gp(mask)).^2,single(0.6)^2.*(f(mask)-g(mask)).^2,'AbsTol',single(1e-4));
verifyEqual(tc,f,single(reshape(1:30,5,6))); % inputs unmodified
end
function testInvalidPathAndCoverage(tc)
mask=true(20); mask(:,10)=false;
verifyFalse(tc,crackPathInDomain([5 5;15 5],mask));
verifyFalse(tc,crackPathInDomain([3 3;15 15;3 15;15 3],true(20)));
verifyFalse(tc,crackPathInDomain([5 5;5 5],true(20)));
R=zeros(20); R(5,5)=NaN;
verifyError(tc,@()phaseFieldEnergy2D(zeros(20),R,true(20),struct()),'phaseFieldEnergy2D:Data');
end
