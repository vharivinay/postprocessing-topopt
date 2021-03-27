%% THIS FUNCTION REDUCES NUMBER OF POINTS OF CURVES FROM BWBOUNDARIES
% Recursive Douglas-Peucker Polyline Simplification

function ps = rdp(p,eps)
%% DEFINE 
last_index     = size(p,1);
first_index     = 1;
% logical vector for the vertices to be retained
I = true(last_index,1);
% call recursive function
p = reduce_line(p,eps,first_index,last_index);
ps = p(I,:);

function p = reduce_line(p,eps,first_index,last_index)
    
    if p(first_index,:)==p(last_index,:)
        % calculate the shortest distance of all vertices
        
        distance = hypot(p(first_index,1)-p(first_index+1:last_index-1,1),...
            p(first_index,2)-p(first_index+1:last_index-1,2));
    else    
        % calculate shortest distance of all points to the line from ixs to ixe
        % subtract starting point from other locations
        pt = p(first_index+1:last_index,:)-p(first_index,:);
        % end point
        a = pt(end,:)';
        beta = (a' * pt')./(a'*a);
        b    = pt-(beta.*a)';
        distance = hypot(b(:,1),b(:,2));
    end
    
    % identify maximum distance and get the linear index of its location
    [dmax,break_index] = max(distance);
    break_index  = first_index + break_index; 
    
    % if the maximum distance is smaller than the tolerance remove vertices
    if dmax <= eps
        if first_index ~= last_index-1
            I(first_index+1:last_index-1) = false;
        end
    % if not, call simplifyrec for the segments between ixs and ixc (ixc
    % and ixe)
    else   
        p = reduce_line(p,eps,first_index,break_index);
        p = reduce_line(p,eps,break_index,last_index);
    end
end

end