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
end
