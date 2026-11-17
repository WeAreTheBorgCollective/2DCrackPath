function cases = screeningCases
% Edit image lists and mesh definitions here; the first image is the reference.
images = arrayfun(@(k) sprintf('Slike/composite%02d.bmp',k), ...
    0:11,'UniformOutput',false);
cases = struct;
cases.original = struct('name','original', ...
    'images',{{'Slike/reference.tif','Slike/deformed.tif'}}, ...
    'meshFile','Mesh/CrackMesh1.m','meshOffset','original');
for suffix = {'050','100','200'}
    name = ['composite' suffix{1}];
    cases.(name) = struct('name',name,'images',{images}, ...
        'meshFile',['Mesh/mesh' suffix{1} '.mat'],'meshOffset',[0 0]);
end
% Windmill experiment: actual acquisition IDs, in temporal order.
windmillIDs = [1 50:50:950 960];
windmillImages = arrayfun(@(k) sprintf('Windmill/windimill_meso_spec_1_%03d.tif',k), ...
    windmillIDs,'UniformOutput',false);
cases.windmill = struct('name','windmill','images',{windmillImages}, ...
    'meshFile','Mesh/Mesh.mat','meshOffset',[0 0]);
end
