function [deformedXY,displacementXY] = plotCrackOnDeformed(ax,g,referenceXY,Mesh,a,label)
%PLOTCRACKONDEFORMED Show the mapped crack on the ORIGINAL graylevel g image.
% No damage and no image warp are applied to the displayed background.
[deformedXY,displacementXY] = mapCrackPathToDeformed(referenceXY,Mesh,a);
imagesc(ax,g);
colormap(ax,gray(256));
axis(ax,'image');
hold(ax,'on');
plot(ax,deformedXY(:,1),deformedXY(:,2),'r-','LineWidth',1.5);
bounds = [min(deformedXY,[],1)-20; max(deformedXY,[],1)+20];
xlim(ax,[max(1,bounds(1,1)) min(size(g,2),bounds(2,1))]);
ylim(ax,[max(1,bounds(1,2)) min(size(g,1),bounds(2,2))]);
title(ax,label,'Interpreter','none');
xlabel(ax,'Deformed x [px]');
ylabel(ax,'Deformed y [px]');
end
