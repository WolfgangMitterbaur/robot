% function to create a cureved linear interpolation between b and c
% Wolfgang Mitterbaur

%% function Interpol2
% input values:
% b: the first point
% c: the second point
% e: the already existing waypoint
% no_points: the number of points for the interpolation
% return value:
% newpoint: all interpolation points

function newpoint = Interpol2(b, c, e, no_points)

    % linear interpolation
    % intermediate points

    t = 1:1:no_points; 
    newpoint = zeros(3, no_points);  % new(dimension, anzahl)
    
    delta_y = c(2) - b(2);  % minimal distance in y
    delta_z = c(3) - b(3);  % minimal distance in z
    delta_y1 = delta_y / (t(end)+1);    % delta distance in y
    delta_z1 = delta_z / (t(end)+1);    % delta distance in z
    
    % start condition
    newpoint(:, 1) = b;
    
    for i = 1:t(end)
        newpoint(1, i) = b(1);
        if i == 1
            newpoint(2,i) = newpoint(2,1) + delta_y1;  % y-value increase
            newpoint(3,i) = newpoint(3,1) + delta_z1;  % z-value increase
        else
            newpoint(2,i) = newpoint(2,i-1) + delta_y1;  % y-value increase
            newpoint(3,i) = newpoint(3,i-1) + delta_z1;  % z-value increase
        end
    end

    delta_offset_y = (newpoint(2,3) - e(2)) / 5;
    delta_offset_z = (newpoint(3,3) - e(3)) / 4;
     
    offset_z = 0;
    damping = 1.0;
    % add a curvature
    for i = 1:t(end)
        if i <= t(end)/2 + 1
            damping = damping - 0.300;
            offset_z = offset_z - delta_offset_z * damping;
        else
            damping = damping + 0.300;
            offset_z = offset_z + delta_offset_z * damping;
        end
        newpoint(3,i) = newpoint(3,i) + offset_z;
    end   
    
    offset_y = 0;
    damping = 1.0;
    % add a curvature
    for i = 1:t(end)
        if i <= t(end)/2 + 2
            damping = damping - 0.30;
            offset_y = offset_y - delta_offset_y  * damping;
        else
            damping = damping + 0.30;
            offset_y = offset_y + delta_offset_y  * damping;
        end
        newpoint(2,i) = newpoint(2,i) + offset_y;
    end 
end
