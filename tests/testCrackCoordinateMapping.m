function tests = testCrackCoordinateMapping
tests = functiontests(localfunctions);
end

function setupOnce(~)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'lib','mesh'));
end

function Mesh = squareMesh
Mesh.pos = [1 1 0;11 1 0;11 11 0;1 11 0];
Mesh.con = [0 1 2 3;0 1 3 4];
end

function testTranslationSignAndNoMutation(tc)
Mesh = squareMesh;
path = [2 3;5 5;8 9];
a = repmat([2 -1 0],4,1);
[mapped,u] = mapCrackPathToDeformed(path,Mesh,a);
verifyEqual(tc,mapped,path+[2 -1],'AbsTol',1e-12);
verifyEqual(tc,u,repmat([2 -1],3,1),'AbsTol',1e-12);
verifyEqual(tc,path,[2 3;5 5;8 9]);
end

function testAffineDisplacement(tc)
Mesh = squareMesh;
B = [0.1 -0.2;0.3 0.05];
translation = [2 -1];
a = [Mesh.pos(:,1:2)*B+translation zeros(4,1)];
path = [2 3;5 5;8 9];
[mapped,u] = mapCrackPathToDeformed(path,Mesh,a);
verifyEqual(tc,u,path*B+translation,'AbsTol',1e-12);
verifyEqual(tc,mapped,path+path*B+translation,'AbsTol',1e-12);
end

function testActualElementConnectivity(tc)
Mesh = squareMesh;
a = zeros(4,3);
a(3,1) = 4;
[~,u] = mapCrackPathToDeformed([8.5 3.5;3.5 8.5],Mesh,a);
verifyEqual(tc,u,[1 0;1 0],'AbsTol',1e-12);
end

function testMeshGapDoesNotExtrapolate(tc)
Mesh.pos = [1 1 0;3 1 0;1 3 0;7 1 0;9 1 0;9 3 0];
Mesh.con = [0 1 2 3;0 4 5 6];
verifyError(tc,@()mapCrackPathToDeformed([5 2],Mesh,zeros(6,3)), ...
    'mapCrackPathToDeformed:OutsideMesh');
end

function testUnsupportedMeshRejected(tc)
Mesh = squareMesh;
Mesh.con(:,1) = 1;
verifyError(tc,@()mapCrackPathToDeformed([5 5],Mesh,zeros(4,3)), ...
    'mapCrackPathToDeformed:Mesh');
end
