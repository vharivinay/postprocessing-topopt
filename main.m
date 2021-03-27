%% MAIN CODE STRUCTURAL DESIGN OPTIMIZATION
tic;
set(0,'DefaultAxesColor','none');

%% ADD PATH TO FUNCTIONS
addpath('/home/harivinay/Desktop/Results/report_code/thesis_func');
%% CODE START -- INPUTS
% Type of Specimen: Allowed specimens are
% =======================================================================
% 'MBB','half-wheel','short_cantilever_endload','short_cantilever_midload',
% 'multiload_cantilever','plate_w_hole_endload','plate_w_hole_midload',
% 'l-beam','t-beam','t-beam_w_hole'
% =======================================================================
specimen = 'MBB';
ft=3; % Filtering Method 1 = sensitivity,2 = desnity or 3 = heaviside
sf = 1.0; % Scaling Factor for dimensions [RECOMMENDED: INTEGERS]
%========================================================================
% validate = 'Y' runs validation routine -  computes compliance & stress
% validate = 'N' Skips validation routine
validate = 'Y';
%========================================================================
% Function call for getting input values
fprintf('Getting input data...\n');
[nelx,nely,volfrac,penal,rmin] = inputdata(specimen,sf);
nG = 8/sf;
fprintf('Done!\n');

%% Calling Function and Importing Data % dens = xPhys, % compliance = c
fprintf('Phase I: Topology Optimization\n')
[xPhys,compliance_unsmooth,...
    beta,vMS] = top88(nelx,nely,volfrac,penal,rmin,ft,specimen);
dens = xPhys;
fprintf('Done!\n');

%% Calling Smooth Function
fprintf('Phase II: Geometry Smoothing\n');
[sdens] = smoothing(dens,nelx,nely,rmin,nG,beta);
% irho = INTERPOLATED and THRESHOLDED DENSITY
fprintf('Phase II: Validation (Optional)\n');
if validate == 'Y'
    [c_smooth,vMS_smooth,svolfrac]=validation(sdens,specimen,nG);
end
fprintf('Done!\n');

%% Boundry Extraction & Smoothing
I = sdens;
fprintf('Phase III: Boundary Extraction\n')
po = 4; % Polynomial Order
[cellsize,smooth_curves,~,~] = extract_boundary(I,po,nelx,nely,nG);
%% Write Smooth Co-ordinate data to File
writecurve(cellsize,smooth_curves);
toc;

%% IN CASE OF BAD RESULTS DUE TO LOW VOLUME FRACTION OR HIGH GRAY VALUES
% Evaluate the following line seperately
% [cu_new,vfu_new,cs_new,vfs_new] = extreme_case(xPhys,beta,specimen,nG);