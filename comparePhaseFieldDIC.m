%% Compare initial and final PF displacements on ORIGINAL unmasked images
% Run after runPhaseFieldDIC. Does not run correlation again or use Damage.
projectRoot=fileparts(mfilename('fullpath')); addpath(projectRoot); setupPaths;
if exist('correli_pathdef','file')==2, addpath(correli_pathdef); end
initial=load(fullfile(projectRoot,'InitialValues','DisplacementInit.mat'));
final=load(fullfile(projectRoot,'Damage','PhaseFieldResult.mat'));
if ~isequal(initial.Mesh.pos,final.Mesh.pos) || ...
        ~isequal(initial.Mesh.con,final.Mesh.con) || ~isequal(size(initial.a),size(final.a))
    error('comparePhaseFieldDIC:Mesh','Initial and final displacement meshes must match.');
end
f=single(imread(fullfile(projectRoot,'Slike','reference.tif')));
g=single(imread(fullfile(projectRoot,'Slike','deformed.tif')));
if ~isequal(size(f),size(g),initial.imageSize)
    error('comparePhaseFieldDIC:ImageSize','Image sizes must match initialization.');
end
Mesh=final.Mesh;
kernel=kernel_corr('cubic',kernel_prep(Mesh,size(f)));
% Native Correli residual = f(X)-g(X+u(X)), using ORIGINAL f and g.
% Thus f-R is the deformed image sampled back on the reference grid.
RInitial=kernel(f,g,reshape(initial.a.',[],1),'R');
RFinal=kernel(f,g,reshape(final.a.',[],1),'R');
warpedInitial=double(f)-double(RInitial);
warpedFinal=double(f)-double(RFinal);
measurementMask=initial.measurementMask;
comparisonMask=measurementMask & isfinite(RInitial) & isfinite(RFinal);
if ~any(comparisonMask(:)), error('comparePhaseFieldDIC:ROI','No common valid pixels.'); end
coverage=[nnz(measurementMask & isfinite(RInitial));nnz(measurementMask & isfinite(RFinal))];
rmse=[sqrt(mean(double(RInitial(comparisonMask)).^2)); ...
      sqrt(mean(double(RFinal(comparisonMask)).^2))];
mae=zeros(2,1);
mae(1)=mean(abs(double(RInitial(comparisonMask))));
mae(2)=mean(abs(double(RFinal(comparisonMask))));
metrics=table(["Initial DIC";"Final PF DIC"],rmse,mae,coverage, ...
    'VariableNames',{'Solution','RMSE_GL','MAE_GL','ValidROIPixels'});
disp(metrics);
deltaA=final.a-initial.a;
fprintf('Common ROI: %d / %d pixels; max nodal displacement change: %.9g px\n', ...
    nnz(comparisonMask),nnz(measurementMask),max(vecnorm(deltaA(:,1:2),2,2)));
if nnz(comparisonMask)<nnz(measurementMask)
    warning('comparePhaseFieldDIC:Coverage','Metrics use the common valid ROI; inspect coverage counts.');
end
% Persist the FINAL ACCEPTED nodal displacement, with compatible init fields.
a=final.a; DisplacementInitPF=a; R=RFinal; imageSize=size(f); grayScale=initial.grayScale;
initialFolder=fullfile(projectRoot,'InitialValues');
save(fullfile(initialFolder,'DisplacementInitPF.mat'),'a','DisplacementInitPF','Mesh', ...
    'R','measurementMask','comparisonMask','imageSize','grayScale','metrics','-v7');
outputFolder=fullfile(projectRoot,'Comparisons');
if ~isfolder(outputFolder), mkdir(outputFolder); end
save(fullfile(outputFolder,'PhaseFieldComparison.mat'),'RInitial','RFinal', ...
    'warpedInitial','warpedFinal','comparisonMask','metrics','deltaA','-v7');

%% Original-image warps and signed residuals with shared scales
fig=figure('Color','w','Position',[40 40 1500 850]);
tiledlayout(fig,2,3,'TileSpacing','compact');
grayLimits=double([min(f(comparisonMask)) max(f(comparisonMask))]);
if grayLimits(1)==grayLimits(2), grayLimits=grayLimits+[-1 1]; end
residualLimit=max([abs(double(RInitial(comparisonMask)));abs(double(RFinal(comparisonMask)));eps]);
showImage(nexttile,double(f),comparisonMask,grayLimits,'Original reference [GL]',true);
showImage(nexttile,warpedInitial,comparisonMask,grayLimits,'Deformed warped by initial u [GL]',true);
showImage(nexttile,warpedFinal,comparisonMask,grayLimits,'Deformed warped by final PF u [GL]',true);
showImage(nexttile,RInitial,comparisonMask,[-residualLimit residualLimit], ...
    sprintf('Initial residual; RMSE %.3f GL',rmse(1)),false);
showImage(nexttile,RFinal,comparisonMask,[-residualLimit residualLimit], ...
    sprintf('Final PF residual; RMSE %.3f GL',rmse(2)),false);
deltaR=double(RFinal)-double(RInitial);
limit=max(abs(deltaR(comparisonMask))); if limit==0, limit=1; end
showImage(nexttile,deltaR,comparisonMask,[-limit limit],'Residual difference: PF - initial [GL]',false);
exportgraphics(fig,fullfile(outputFolder,'ResidualComparison.png'),'Resolution',160);

%% Displacements on the actual finite element mesh, shared component scales
fig=figure('Color','w','Position',[40 40 1500 850]);
tiledlayout(fig,2,3,'TileSpacing','compact');
for component=1:2
    values=[initial.a(:,component);final.a(:,component)];
    limits=[min(values) max(values)]; if limits(1)==limits(2), limits=limits+[-1 1]; end
    label='x'; if component==2, label='y'; end
    showDisplacement(nexttile,Mesh,initial.a(:,component),limits,['Initial u_' label ' [px]']);
    showDisplacement(nexttile,Mesh,final.a(:,component),limits,['Final PF u_' label ' [px]']);
    bound=max(abs(deltaA(:,component))); if bound==0, bound=1; end
    showDisplacement(nexttile,Mesh,deltaA(:,component),[-bound bound],['Change in u_' label ' [px]']);
end
exportgraphics(fig,fullfile(outputFolder,'DisplacementComparison.png'),'Resolution',160);

function showImage(ax,image,mask,limits,label,useGray)
im=imagesc(ax,image,limits); im.AlphaData=mask;
axis(ax,'image'); colorbar(ax); title(ax,label,'Interpreter','none','FontSize',13);
set(ax,'FontSize',10); xlabel(ax,'x [px]'); ylabel(ax,'y [px]');
if useGray, colormap(ax,gray); else, colormap(ax,turbo); end
end
function showDisplacement(ax,Mesh,values,limits,label)
patch(ax,'Faces',Mesh.con(:,2:4),'Vertices',Mesh.pos(:,1:2), ...
    'FaceVertexCData',values,'FaceColor','interp','EdgeColor','none');
axis(ax,'image'); set(ax,'YDir','reverse','FontSize',10); clim(ax,limits);
colorbar(ax); colormap(ax,turbo); title(ax,label,'Interpreter','none','FontSize',13);
xlabel(ax,'x [px]'); ylabel(ax,'y [px]');
end
