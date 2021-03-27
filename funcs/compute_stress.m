function vMS = compute_stress(nelx,nely,penal,specimen,F,U,edofMat,dens)
E0 = 1;
Emin = 1e-9;
nu = 0.3;
L = 1;
%% PREPARE STRESS ANALYSIS
B = (1/2/L)*[-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1; -1 -1 -1 1 1 1 1 -1];
DE = (1/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]; 
E = Emin+dens(:)'.^penal*(E0-Emin);
% STRESS FOR MULTILOAD CASE
ss = zeros(3*nelx*nely,1);
s1 = reshape(ss,nelx*nely,3);
s2 = reshape(ss,nelx*nely,3);

%% STRESS CALCULATION
switch specimen
    case {'multiload_cantilever','t-beam','t-beam_w_hole','c-bracket'}
        for i = 1:size(F,2)
            if i == 1
              Ui = U(:,i);
              s1 = (Ui(edofMat)*(DE*B)').*repmat(E',1,3);
            elseif i == 2 
              Ui = U(:,i);
              s2 = (Ui(edofMat)*(DE*B)').*repmat(E',1,3);
            end
        end
        s = s1+s2;
        vMS = reshape(sqrt(s(:,1).^2+s(:,2).^2-s(:,1).*s(:,2)+3.*s(:,3).^2),nely,nelx);
    otherwise
    s = (U(edofMat)*(DE*B)').*repmat(E',1,3);
        vMS = reshape(sqrt(s(:,1).^2+s(:,2).^2-s(:,1).*s(:,2)+3.*s(:,3).^2),nely,nelx);

end