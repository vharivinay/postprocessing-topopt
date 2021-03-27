%% THIS FUNCTION RESIZES RHO TO AVOID DECIMELS WHILE APPLYTING BOUNDARY CONDITIONS
function sdens = resize_rho(rho)
    I = mat2gray(rho);
    [y,x] = size(I);
    I = imresize(I,[y-1 x-1]);
    [p,q] = size(I);
    sdens = zeros(p,q);
    for i = 1:p
        for j = 1:q
           sdens(i,j) = I(i,j);
        end
    end
    sdens = (sdens - min(sdens))./(max(sdens)-min(sdens)); 
end