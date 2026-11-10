function tests=testResidualCrackPath
tests=functiontests(localfunctions);
end
function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'lib','optimization'));
end
function testNgEquationAndBinBoundary(tc)
% Direct evaluation of Eq.8, independently forming the two classes.
counts=[30 40 20 2 1 1 5 8];
v=repelem(0:7,counts); p=counts/sum(counts);
score=-inf(1,8);
for t=0:6
    a=0:t; b=t+1:7;
    wa=sum(p(a+1)); wb=sum(p(b+1));
    ma=sum(a.*p(a+1))/wa; mb=sum(b.*p(b+1))/wb;
    score(t+1)=(1-p(t+1))*(wa*ma^2+wb*mb^2);
end
[~,k]=max(score);
[bw,T,info]=valleyOtsuThreshold([v NaN Inf],8);
verifyEqual(tc,T,k-0.5,'AbsTol',1e-12);
verifyEqual(tc,info.objective(1:7),score(1:7)','AbsTol',1e-12);
verifyEqual(tc,bw,[v>=T false false]);
% Multiplication preserves zero-based levels and scales the GL threshold.
[scaled,Ts]=valleyOtsuThreshold(10*v,8);
verifyEqual(tc,scaled,bw(1:end-2)); verifyEqual(tc,Ts,10*T,'AbsTol',1e-12);
end
function testThresholdDegenerate(tc)
[bw,T]=valleyOtsuThreshold(nan(5)); verifyFalse(tc,any(bw(:))); verifyTrue(tc,isnan(T));
[bw,T]=valleyOtsuThreshold(7*ones(5)); verifyFalse(tc,any(bw(:))); verifyEqual(tc,T,7);
[bw,~,info]=valleyOtsuThreshold([zeros(1,100) ones(1,2)],256);
verifyEqual(tc,nnz(bw),2); verifyEqual(tc,sum(info.probability),1,'AbsTol',1e-12);
verifyError(tc,@()valleyOtsuThreshold([-1 2]),'valleyOtsuThreshold:Negative');
end
function testThinningTopology(tc)
bw=false(30,50); bw(10:16,3:45)=true;
s=thinCrackMask(bw);
verifyFalse(tc,any(s(:)&~bw(:)));
cc=bwconncomp(s,8); verifyEqual(tc,cc.NumObjects,1);
verifyFalse(tc,any(conv2(double(s),ones(2),'valid')==4,'all'));
verifyEqual(tc,thinCrackMask(s),s);
verifyFalse(tc,any(thinCrackMask(false(5)),'all'));
end
function testFragmentConnection(tc)
s=false(30,50); s(10,2:10)=true; s(10,15:22)=true; s(10,30:40)=true;
o=struct('MinComponentDiameter',0,'MaxVertexDensity',Inf,'MaxGap',10);
[xy,mask,G,info]=connectCrackGraph(s,o);
verifyEqual(tc,info.bridgeCount,2); verifyTrue(tc,info.allRetainedConnected);
verifyEqual(tc,info.pathLength,38,'AbsTol',1e-12);
verifyEqual(tc,xy([1 end],:),[2 10;40 10]);
verifyTrue(tc,all(mask(10,2:40)));
verifyEqual(tc,numel(unique(conncomp(G))),1);
verifyEqual(tc,G.Edges.Weight(G.Edges.IsBridge),[5;8],'AbsTol',1e-12);
end
function testGapLimitAndInvalidRegion(tc)
s=false(20,30); s(10,2:10)=true; s(10,15:25)=true;
o=struct('MinComponentDiameter',0,'MaxVertexDensity',Inf,'MaxGap',4);
[~,~,~,info]=connectCrackGraph(s,o);
verifyEqual(tc,info.remainingComponents,2); verifyEqual(tc,info.pathLength,10);
o.MaxGap=100; o.ValidMask=true(size(s)); o.ValidMask(:,12)=false;
[~,mask,~,info]=connectCrackGraph(s,o);
verifyEqual(tc,info.bridgeCount,0); verifyFalse(tc,any(mask(:,12)));
end
function testEuclideanWeightsAndBranchChoice(tc)
s=false(20); s(10,2:18)=true; s(6:9,10)=true;
o=struct('MinComponentDiameter',0,'MaxVertexDensity',Inf,'MaxGap',0);
[xy,~,G,info]=connectCrackGraph(s,o);
verifyEqual(tc,xy([1 end],:),[2 10;18 10]);
verifyEqual(tc,info.pathLength,16,'AbsTol',1e-12);
ends=G.Edges.EndNodes;
expected=hypot(G.Nodes.X(ends(:,1))-G.Nodes.X(ends(:,2)), ...
    G.Nodes.Y(ends(:,1))-G.Nodes.Y(ends(:,2)));
verifyEqual(tc,G.Edges.Weight,expected,'AbsTol',1e-12);
verifyGreaterThan(tc,info.excludedNodeCount,0);
verifyEqual(tc,numel(unique(info.pathNodeIDs)),numel(info.pathNodeIDs));
end
function testEmptySingletonLoopAndFiltering(tc)
[xy,mask,G]=connectCrackGraph(false(10));
verifyEmpty(tc,xy); verifyFalse(tc,any(mask(:))); verifyEqual(tc,numnodes(G),0);
s=false(10); s(5,5)=true;
o=struct('MinComponentDiameter',0,'MaxVertexDensity',Inf,'MaxGap',0);
[xy,mask,~,info]=connectCrackGraph(s,o);
verifyEqual(tc,xy,[5 5]); verifyEqual(tc,nnz(mask),1); verifyEqual(tc,info.pathLength,0);
s=false(10); s(3,3:7)=true; s(7,3:7)=true; s(3:7,3)=true; s(3:7,7)=true;
[xy,~,~,info]=connectCrackGraph(s,o);
verifyGreaterThan(tc,size(xy,1),1);
verifyEqual(tc,numel(unique(info.pathNodeIDs)),numel(info.pathNodeIDs));
[xy,~,G]=connectCrackGraph(s); verifyEmpty(tc,xy); verifyEqual(tc,numnodes(G),0);
end
