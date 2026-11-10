function tests=testRethoreLevelSets
tests=functiontests(localfunctions);
end
function setupOnce(~)
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib','optimization'));
end
function testPlanarPropagationBothDirections(tc)
[~,y]=meshgrid(1:100,1:100); phi=y-50;
for speed=[-3 3]
 actual=advectLevelSetGodunov(phi,speed*ones(size(phi)),1);
 verifyEqual(tc,actual(40:60,40:60),phi(40:60,40:60)-speed,'AbsTol',1e-10);
end
verifyEqual(tc,advectLevelSetGodunov(phi,zeros(size(phi)),1),phi);
end
function testNegativeResidualRidgeAndTipRebuild(tc)
[~,y]=meshgrid(1:60,1:50); R=-100*exp(-(y-24).^2/2);
LS=crackLevelSets2D([(10:50)' 20*ones(41,1)],size(R));
[new,d,info]=updateCrackLevelSetsRethore(LS,R,true(size(R)),struct());
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
verifyFalse(tc,info.accepted); verifyEqual(tc,new,LS);
end
function testCannotJumpDomainHole(tc)
[~,y]=meshgrid(1:60,1:50); R=100*exp(-(y-27).^2/2);
mask=true(size(R)); mask(23:24,:)=false; R(~mask)=NaN;
LS=crackLevelSets2D([(10:50)' 20*ones(41,1)],size(R));
[new,~,~]=updateCrackLevelSetsRethore(LS,R,mask,struct());
verifyLessThanOrEqual(tc,max(new.pathXY(:,2)),22+1e-8);
verifyTrue(tc,crackPathInDomain(new.pathXY,mask));
end
