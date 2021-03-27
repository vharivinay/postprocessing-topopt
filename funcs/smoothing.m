function [sdens] = smoothing(xPhys,nelx,nely,rmin,nG,beta)

%% CREATE FINER (INTERPOLATED) MESHGRID
[nodex,nodey] = meshgrid(0:nelx,0:nely);
[fnx,fny] = meshgrid(0:1/nG:nelx,0:1/nG:nely);
%% PREPARE FILTER FOR NODAL DENSITIES
inH = ones((nelx+1)*(nely+1)*(2*(ceil(rmin)+1))^2,1) ;
jnH = ones(size(inH)); 
snH = zeros(size(inH));
k =0;
for in1 =1:nelx+1
    for jn1 = 1:nely+1
        en1 = (in1-1)*(nely+1)+jn1;
        for in2 = max(in1-ceil(rmin),1):min(in1+ceil(rmin)-1,nelx)
            for jn2 = max(jn1-ceil(rmin),1):min(jn1+ceil(rmin)-1,nely)
                en2 = (in2-1)*nely+jn2;
                k = k+1;
                inH(k) = en1 ;
                jnH(k) = en2 ;
                snH(k) = max(0,rmin-sqrt((in1-in2)^2+(jn1-jn2)^2));
            end
        end
    end
end

Hn = sparse(inH,jnH,snH);
Hns = sum(Hn,2);
%% ASSIGN FILTERED ELEMENTAL VOLUME FRACTIONS TO NODAL DENSITIES
ndens = reshape((Hn*xPhys(:)./Hns),nely+1,nelx+1);

%% UPDATE NODAL DESNITY BY A HEAVISIDE SMOOTH FUNCTION
% fprintf('Running Interpolation...\n')
idens = interp2(nodex,nodey,ndens,fnx,fny,'linear');
eta = (max(idens(:))+min(idens(:)))/2;
irho = max(0,(tanh(beta*eta)+tanh(beta*(idens-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta))));

%% CONVERT TO GRAY SCALE IMAGE
sdens = mat2gray(irho);

%% RESIZE RHO TO AVOID FLOATING POINT VALUES FOR BCs
% fprintf('Preparing density values...\n')
[ys,xs] = size(sdens);
sdens = imresize(sdens,[ys-1 xs-1]);

end