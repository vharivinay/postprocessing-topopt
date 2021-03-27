function [cu_new,vfu_new,cs_new,vfs_new] = extreme_case(dens,beta,specimen,nG)
[nelx,nely,~,~,rmin] = inputdata(specimen);

eta = (min(dens(:))+max(dens(:)))/2;
dens = max(0,(tanh(beta*eta)+tanh(beta*(dens-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta))));

[cu_new,~,vfu_new]=validation(dens,specimen,1);
%% Calling Smooth Function
fprintf('Phase II: Geometry Smoothing\n');
[sdens] = smoothing(dens,nelx,nely,rmin,nG,beta);
% irho = INTERPOLATED and THRESHOLDED DENSITY
[cs_new,~,vfs_new]=validation(sdens,specimen,nG);
fprintf('Done!\n');

thval = graythresh(sdens);
sdens = imbinarize(sdens,thval);

%% Boundry Extraction & Smoothing
I = sdens;
fprintf('Phase III: Boundary Extraction\n')
po = 4; % Polynomial Order
[cellsize,smooth_curves,~,~] = extract_boundary(I,po,nelx,nely,nG);
%% Write Smooth Co-ordinate data to File
writecurve(cellsize,smooth_curves);
end