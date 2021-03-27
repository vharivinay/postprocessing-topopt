%% THIS FUNCTION DEFINES PASSIVE ELEMENTS FOR SPECIMENS WITH FIXED VOIDS
function passive = void_elements(nelx,nely,specimen)
    switch specimen
        case {'plate_w_hole_endload','plate_w_hole_midload'}
            passive = zeros(nely,nelx);
            for i = 1:nelx
                for j = 1:nely
                    if sqrt((j-nely/2)^2+(i-nelx/3)^2) < nely/3
                    passive(j,i) = 1;
                    end
                end
            end
        case 'l-beam'
            passive = zeros(nely,nelx);
            for i = 1:nelx
                for j = 1:nely
                    if i > (0.4*nelx) && j < ((1-0.4)*nely)
                        passive(j,i) = 1;
                    end
                end
            end
        case 't-beam'
            passive = zeros(nely,nelx);
            for i = 1:nelx
                for j = 1:nely
                    if (i <= (0.35*nelx) || i > (0.65*nelx)) && j < ((1-0.3)*nely)
                        passive(j,i) = 1;
                    end
                end
            end
        case 't-beam_w_hole'
            passive = zeros(nely,nelx);
            for i = 1:nelx
                for j = 1:nely
                    if (i <= (0.35*nelx) || i > (0.65*nelx)) && j < ((1-0.3)*nely) ||...
                            i == 0.81*nely && j == 0.85*nely...
                            || sqrt((j-(0.85*nely))^2+(i-(0.75*nelx))^2) < nely/12
                        passive(j,i) = 1;
                    end
                end
            end
        case 'c-bracket'
            passive = zeros(nely,nelx);
            for i = 1:nelx
                for j = 1:nely
                    if i > (0.5*nelx)
                        if j > (0.4*nely)
                            if j < (0.6*nely)
                                passive(j,i) = 1;
                            end
                        end
                    end
                end
            end
    end
end