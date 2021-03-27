function [c,vMS,dc,U] = fe_compliance(xs,ys,penal,specimen,sdens)
nelx = xs;
nely = ys;
E0 = 1;
Emin = 1e-9;
nu = 0.3;
xPhys = sdens;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
[loads,displacements,constraints]=boundary_conditions(nelx,nely,specimen);
F = loads;
U = displacements;
fixeddofs = constraints;
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% FE-ANALYSIS
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
c=0;
dc=0;
switch specimen
  case {'multiload_cantilever','t-beam','t-beam_w_hole','c-bracket'}
      U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
      %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      for i = 1:size(F,2)
          Ui = U(:,i);
          ce = reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
          c = c+ sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
          dc=dc-penal*(E0-Emin)*xPhys.^(penal-1).*ce;
      end
  otherwise
      U(freedofs) = K(freedofs,freedofs)\F(freedofs);
      %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
      c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
      dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
end
 vMS = compute_stress(nelx,nely,penal,specimen,F,U,edofMat,xPhys);
end