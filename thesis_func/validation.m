function [c_smooth,vMS_smooth,fvolfrac]=validation(sdens,specimen,nG)
%% Get necessary Inputs
[~,~,~,penal,~] = inputdata(specimen);

%% CHECK FOR UNDEFINED VALUES AND REPLACE THEM APPROPRIATELY WITH 0/1 
sdens(isnan(sdens)) = mean(sdens(:));
sdens(isinf(sdens)) = mean(sdens(:));
sdens = real(sdens);

[ys,xs]=size(sdens);
%% COMPUTE COMPLIANCE & STRESS FOR SMOOTH GEOMETRY
[c_smooth,vMS_smooth,~] = fe_compliance(xs,ys,penal,specimen,sdens);

%% COMPUTE FINAL VOLFRAC FOR SMOOTH GEOMETRY
fvolfrac = mean(sdens(:));

%% VISUALIZE RESULTS
figure(2);
ax2 = subplot(1,2,1);
cMap = gray;
imagesc((1-sdens)); caxis([0 1]);
caption1 = sprintf('Smooth geometry \n Volfrac = %d',fvolfrac);
title(caption1);
colormap(ax2,cMap);axis equal; axis off; drawnow;
hold on;

%% ADDING COMPLIANCE VALUE IN THE PLOT
caption1 = sprintf('Smooth geometry \n Volfrac = %d \n Compliance = %d',fvolfrac,c_smooth);
title(caption1);
hold off;

%% STRESS DISTRIBUTION OF SMOOTH GEOMETRY
ax1 = subplot(1,2,2);
colormap;imagesc(nG*vMS_smooth); caxis([0 1]);
caption2 = sprintf('von Mises Stress distribution');
title({caption2,' '});
axis equal off;
cMap = jet(256);
cMap(1,:) = 1;
colormap(ax1,cMap);colorbar;drawnow;
end