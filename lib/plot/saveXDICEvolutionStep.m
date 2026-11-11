function saveXDICEvolutionStep(folder,step,LS,solution,R,prep,info,g,Mesh)
%SAVEXDICEVOLUTIONSTEP Persist all three level sets and both physical crack faces.
[minusFace,plusFace] = mapXDICCrackFaces(LS.pathXY,Mesh,solution,prep);
filename = fullfile(folder,sprintf('Step_%02d',step));
save([filename '.mat'],'LS','solution','R','prep','info','minusFace','plusFace','-v7.3');
fig = figure('Visible','off','Color','w','Position',[40 40 1400 850]);
cleanup = onCleanup(@()close(fig)); %#ok<NASGU>
bounds = [min(LS.pathXY,[],1)-15;max(LS.pathXY,[],1)+15];
tiledlayout(fig,3,1,'TileSpacing','compact');
ax = nexttile;
h = imagesc(ax,abs(R));
h.AlphaData = isfinite(R);
hold(ax,'on');
plot(ax,LS.pathXY(:,1),LS.pathXY(:,2),'r-');
axis(ax,'image');
xlim(ax,bounds(:,1)');
ylim(ax,bounds(:,2)');
colorbar(ax);
title(ax,sprintf('Step %d: original-image |f(X)-g(X+u)|',step));
ax = nexttile;
imagesc(ax,LS.phi);
hold(ax,'on');
contour(ax,LS.phi,[0 0],'r','LineWidth',1.5);
contour(ax,LS.psiStart,[0 0],'k--','LineWidth',1.5);
contour(ax,LS.psiEnd,[0 0],'k--','LineWidth',1.5);
axis(ax,'image');
xlim(ax,bounds(:,1)');
ylim(ax,bounds(:,2)');
colorbar(ax);
title(ax,'Surface phi=0 (red); endpoint fronts psiStart=0 and psiEnd=0 (dashed)');
ax = nexttile;
imagesc(ax,g);
colormap(ax,gray);
axis(ax,'image');
hold(ax,'on');
plot(ax,minusFace(:,1),minusFace(:,2),'c-','LineWidth',1.5);
plot(ax,plusFace(:,1),plusFace(:,2),'m-','LineWidth',1.5);
center = (minusFace+plusFace)/2;
xlim(ax,[min(center(:,1))-15 max(center(:,1))+15]);
ylim(ax,[min(center(:,2))-15 max(center(:,2))+15]);
title(ax,'Original deformed image: two XFEM crack-face traces (cyan / magenta)');
set(findall(fig,'Type','Axes'),'Toolbar',[]);
exportgraphics(fig,[filename '.png'],'Resolution',140);
end
