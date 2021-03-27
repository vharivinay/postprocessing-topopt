function [cu_new,vfu_new,cs_new,vfs_new] = extreme_case(dens,beta,specimen,nG)
[nelx,nely,volfrac,~,rmin] = inputdata(specimen);

% set(groot,'defaultFigureVisible','off')
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
% move = 0.001;
% ch = 1;
% loop=0;
% 
% while ch > 0.001 && loop < 1000
%     vf_one = volfrac;
%     th_dens = imbinarize(sdens,thval);
%     vf_two = mean(th_dens(:));
%     if vf_two == vf_one
%         break;
%     elseif vf_two > vf_one
%         thval = thval + move;
%     else
%         thval = thval - move;
%     end
%     loop = loop+1;
%     ch = abs(vf_one-vf_two);
% end
% sdens = th_dens;

%% Boundry Extraction & Smoothing
% set(groot,'defaultFigureVisible','on')
I = sdens;
fprintf('Phase III: Boundary Extraction\n')
po = 4; % Polynomial Order
[cellsize,smooth_curves,~,~] = extract_boundary(I,po,nelx,nely,nG);
%% Write Smooth Co-ordinate data to File
writecurve(cellsize,smooth_curves);
end