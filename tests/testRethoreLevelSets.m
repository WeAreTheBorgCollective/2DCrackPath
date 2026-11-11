function tests=testRethoreLevelSets
tests=functiontests(localfunctions);
end
function setupOnce(~)
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib','optimization'));
end
function testPlanarPropagationBothDirections(tc)
[~,y]=meshgrid(1:100,1:100);
phi=y-50;
for speed=[-3 3]
 actual=advectLevelSetGodunov(phi,speed*ones(size(phi)),1);
 verifyEqual(tc,actual(40:60,40:60),phi(40:60,40:60)-speed,'AbsTol',1e-10);
end
verifyEqual(tc,advectLevelSetGodunov(phi,zeros(size(phi)),1),phi);
end
function testNegativeResidualRidgeAndTipRebuild(tc)
[~,y]=meshgrid(1:60,1:50);
R=-100*exp(-(y-24).^2/2);
LS=crackLevelSets2D([(10:50)' 20*ones(41,1)],size(R));
[new,d,info]=updateCrackLevelSetsRethore(LS,R,true(size(R)),struct(),struct('screeningWindowSize',[1 1],'velocityMethod','maximum'));
verifyTrue(tc,info.accepted);
verifyGreaterThan(tc,info.meanResidualAfter,info.meanResidualBefore);
verifyEqual(tc,median(new.pathXY(:,2)),24,'AbsTol',0.3);
verifyTrue(tc,crackPathInDomain(new.pathXY,true(size(R))));
verifyEqual(tc,d,createDamageField(new,true(size(R)),struct()),'AbsTol',1e-12);
rebuilt=crackLevelSets2D(new.pathXY,size(R));
verifyEqual(tc,new.psiStart,rebuilt.psiStart);
verifyEqual(tc,new.psiEnd,rebuilt.psiEnd);
verifyEqual(tc,interp2(new.psiStart,new.pathXY(1,1),new.pathXY(1,2)),0,'AbsTol',1e-10);
verifyEqual(tc,interp2(new.psiEnd,new.pathXY(end,1),new.pathXY(end,2)),0,'AbsTol',1e-10);
end
function testFlatResidualNoDrift(tc)
LS=crackLevelSets2D([10 20;50 20],[50 60]);
[new,~,info]=updateCrackLevelSetsRethore(LS,ones(50,60),true(50,60),struct());
verifyFalse(tc,info.accepted);
verifyEqual(tc,new,LS);
end
function testCannotJumpDomainHole(tc)
[~,y]=meshgrid(1:60,1:50);
R=100*exp(-(y-27).^2/2);
mask=true(size(R));
mask(23:24,:)=false;
R(~mask)=NaN;
LS=crackLevelSets2D([(10:50)' 20*ones(41,1)],size(R));
[new,~,~]=updateCrackLevelSetsRethore(LS,R,mask,struct(),struct('screeningWindowSize',[1 1],'velocityMethod','maximum'));
verifyLessThanOrEqual(tc,max(new.pathXY(:,2)),22+1e-8);
verifyTrue(tc,crackPathInDomain(new.pathXY,mask));
end

function testIndependentScreeningWindows(tc)
R=zeros(80,80);
R(38,:)=100;
R(42:48,:)=30;
LS=crackLevelSets2D([(20:60)' 40*ones(41,1)],size(R));
[~,~,raw]=updateCrackLevelSetsRethore(LS,R,true(size(R)),struct(), ...
    struct('screeningWindowSize',[1 1],'velocityMethod','maximum'));
[new,~,screen]=updateCrackLevelSetsRethore(LS,R,true(size(R)),struct(), ...
    struct('screeningWindowSize',[7 3],'velocityMethod','maximum'));
verifyLessThan(tc,median(raw.targetPoints(:,2)),40);
verifyGreaterThan(tc,median(screen.targetPoints(:,2)),40);
verifyTrue(tc,screen.accepted);
verifyGreaterThan(tc,median(new.pathXY(:,2)),40);
verifyEqual(tc,screen.screeningWindowSize,[7 3]);
verifyGreaterThan(tc,screen.meanScreeningAfter,screen.meanScreeningBefore);
end
function testDefaultWindowAndNoCompleteWindows(tc)
R=ones(50,60);
mask=false(size(R));
mask(20:30,:)=true;
LS=crackLevelSets2D([15 25;45 25],size(R));
[new,~,info]=updateCrackLevelSetsRethore(LS,R,mask,struct());
verifyEqual(tc,info.screeningWindowSize,[23 23]);
verifyEqual(tc,info.validScreeningCenters,0);
verifyFalse(tc,info.accepted);
verifyEqual(tc,new,LS);
verifyEqual(tc,info.velocity,zeros(2,1));
end

function testGradientAscentBothSigns(tc)
[~,y]=meshgrid(1:80,1:80);
LS=crackLevelSets2D([(20:60)' 40*ones(41,1)],size(y));
for direction=[-1 1]
    R=100+direction*y;
    options=struct('screeningWindowSize',[5 7], ...
        'gradientGain',0.5,'maxNormalStep',1);
    [new,~,info]=updateCrackLevelSetsRethore(LS,R,true(size(R)),struct(),options);
    verifyTrue(tc,info.accepted);
    verifyEqual(tc,info.velocityMethod,'gradient');
    verifyEqual(tc,median(new.pathXY(:,2)),40+direction*0.5,'AbsTol',0.05);
    verifyGreaterThan(tc,info.meanScreeningAfter,info.meanScreeningBefore);
end
end

function testGradientStepCap(tc)
[~,y]=meshgrid(1:80,1:80);
LS=crackLevelSets2D([(20:60)' 40*ones(41,1)],size(y));
options=struct('screeningWindowSize',[3 3],'maxNormalStep',0.25);
[~,~,info]=updateCrackLevelSetsRethore(LS,100*y,true(size(y)),struct(),options);
verifyLessThanOrEqual(tc,max(abs(info.velocity)),0.25);
verifyTrue(tc,info.accepted);
end
