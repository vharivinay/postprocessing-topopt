%% THIS FUNCTION CONTAINS DEFINED BOUNDARY CONDITIONS FOR ALL THE CASES
function [F,U,fixeddofs] = boundary_conditions(nelx,nely,specimen)
    switch specimen
        case 'MBB'
            F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
        case 'half-wheel'
            F = sparse(((nely+1)*(nelx+1))+nely+1,1,-1,2*(nely+1)*(nelx+1),1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            fixeddofs = union([2*(nely+1)],[2*(nelx+1)*(nely+1)-1:1:2*(nelx+1)*(nely+1)]);
        case 'short_cantilever_endload'
            F = sparse(2*(nely+1)*(nelx+1),1,-1,2*(nely+1)*(nelx+1),1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            fixeddofs = [2:2*nely+1];
        case 'short_cantilever_midload'
            F = sparse((2*(nely+1)*(nelx+1))-(((2*(0.5*nely+1)))),...
                1,-1,2*(nely+1)*(nelx+1),1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            fixeddofs = [2:2*nely+1];
        case 'multiload_cantilever'
            F = sparse([2*(nely+1)*nelx+2,2*(nely+1)*(nelx+1)],[1 2],[1 -1],...
                2*(nely+1)*(nelx+1),2);
            U = zeros(2*(nely+1)*(nelx+1),2);
            fixeddofs = [1:2*nely+1];
        case 'plate_w_hole_endload'
            F = sparse(2*(nely+1)*(nelx+1),1,-1,2*(nely+1)*(nelx+1),1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            fixeddofs = [2:2*nely+1];
        case 'plate_w_hole_midload'
            F = sparse(2*(nely+1)*(nelx+1)-(((2*((0.5*nely)+2)))),1,-1,2*(nely+1)*(nelx+1),1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            fixeddofs = [2:2*nely+1];
        case 'l-beam'
            F = sparse((2*(nely+1)*(nelx+1)-((2*(0.4*nely+1)))),1,-1,2*(nely+1)*(nelx+1),1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            fixeddofs = union([2*(nelx+1)*(0:nely)+1], [2*(nelx+1)*(0:nely)+2]);
        case 't-beam'
            F = sparse([(2*(0.7*nely)),2*(nely+1)*(nelx+1)-((2*(0.3*nely+1)))],...
                [1 2],[-1 -1],2*(nely+1)*(nelx+1),2);
            U = zeros(2*(nely+1)*(nelx+1),2);
            fixeddofs = union([2*(nely+1)*(0:nelx)+1], [2*(nely+1)*(0:nelx)+2]);
        case 't-beam_w_hole'
            F = sparse([(2*(0.7*nely)),2*(nely+1)*(nelx+1)-((2*(0.3*nely+1)))],...
                [1 2],[-1 -1],2*(nely+1)*(nelx+1),2);
            U = zeros(2*(nely+1)*(nelx+1),2);
            fixeddofs = union([2*(nely+1)*(0:nelx)+1], [2*(nely+1)*(0:nelx)+2]);
        case 'c-bracket'
            F = sparse([(2*(nely+1)*(nelx+1)-(0.62*(2*(nely)))),...
                (2*(nely+1)*(nelx+1)-(0.4*(2*(nely))))],[1 2],[1 -1],2*(nely+1)*(nelx+1),2);
            U = zeros(2*(nely+1)*(nelx+1),2);
            fixeddofs = union([(nely+2)],...
                [(2*0.5*(nely+1)*(nelx+1):1:2*0.5*(nely+1)*(nelx+1)+1)]);
     end
end
