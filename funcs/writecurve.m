%% THIS FUNCTION WRITES CURVES TO TXT WHICH CAN BE IMPORTED IN SPACECLAIM
function writecurve(cellsize,smooth_curves)
for m = 1:cellsize
    thiscellcontents = smooth_curves{m};
    q = size(thiscellcontents,1);
    w = ones(q,1); % Creating a row of ones for the txt file
    x = thiscellcontents(:,2);    
    y = thiscellcontents(:,1);
    digits(4);
    vpa(x);
    vpa(y);
    r = [w x y];
    filename = sprintf('%d.txt', m);
    fid = fopen(filename, 'w');
    fprintf(fid, "Polyline = true\n\n");
    fclose(fid);
    fid = fopen(filename, 'a');
    dlmwrite(filename,r,'-append','delimiter','\t');
    fclose(fid);           
end
end
%% ------------- FOR REUSE IN THE FUTURE---------------------------
% SPACECLAIM POINTS IMPORT FORMAT
%     r = [w x y];
%     filename = sprintf('%d.txt', m);
%     fid = fopen(filename, 'w');
%     fprintf(fid, "Polyline = true\n\n");
%     fclose(fid);
%     fid = fopen(filename, 'a');
%     dlmwrite(filename,r,'-append','delimiter','\t');
%     fclose(fid);

% DXF IMPORT FORMAT %---------REQUIRES dxflib FUNCTION------------
% LINK TO DXFLIB: https://de.mathworks.com/matlabcentral/fileexchange/33884-dxflib
%     r = [x y z];
%     filename = sprintf('%d.dxf', m);
%     fid = dxf_open(filename);
%     dxf_polyline(fid,x,y,x);
%     dxf_close(fid);