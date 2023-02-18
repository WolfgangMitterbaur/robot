% function to create a cureved linear interpolation between b and c
% Wolfgang Mitterbaur

%% function Interpol
% input values:
% a: the previous point
% b: the first point
% c: the second point
% d: the following point
% no_points: the number of points for the interpolation
% return value:
% newpoint: all interpolation points

function newpoint = Interpol(a, b, c, d, no_points)

    % distance between the prevoius point and the first point
    deltay_ab = (b(2) - a(2));  % y-distance between a and b
    deltay_bc = (c(2) - b(2));  % y-distance between b and c
    deltay_cd = (d(2) - c(2)); % y-distance between c and d
    
    deltaz_ab = (b(3) - a(3)); % z-distance bewtween a and b
    deltaz_bc = (c(3) - b(3)); % z-distance between b and c
    deltaz_cd = (d(3) - c(3)); % z-distance between c and d
    
    % grandient bewtween the points
    if deltay_ab ~= 0
        steig_ab = deltaz_ab / deltay_ab;
    else
        if deltaz_ab >= 0
            steig_ab = 99999999;
        else 
            steig_ab = -99999999;
        end
    end
    
    if deltay_bc ~= 0
        steig_bc = deltaz_bc / deltay_bc;
    else
        if deltaz_bc >= 0
            steig_bc = 99999999;
        else 
            steig_bc = -99999999;
        end
    end
    
    if deltay_cd ~= 0
        steig_cd = deltaz_cd / deltay_cd;
    else
        if deltaz_cd >= 0
            steig_cd = 99999999;
        else 
            steig_cd = -99999999;
        end
    end
    
    % curve from left to right
    % y-coordinates getting bigger
    left_right_h1 = ((a(2) < b(2) && b(2) <= c(2) && c(2) <= d(2))) || ...
        ((a(2) <= b(2) && b(2) < c(2) && c(2) <= d(2))) || ...
        ((a(2) <= b(2) && b(2) <= c(2) && c(2) < d(2)));
    
    left_right_n1 = ((a(2) < b(2) && b(2) <= c(2))) || ...
        ((a(2) <= b(2) && b(2) < c(2)));
    
    
    left_right_n2 = ((b(2) < c(2) && c(2) <= d(2))) || ...
        ((b(2) <= c(2) && c(2) < d(2)));
    
    % curve from right to left
    % y-coordinates getting smaller
    right_left_h1 = ((a(2) > b(2) && b(2) >= c(2) && c(2) >= d(2))) || ...
        ((a(2) >= b(2) && b(2) > c(2) && c(2) >= d(2))) || ...bottom_top_n2
        ((a(2) >= b(2) && b(2) >= c(2) && c(2) > d(2)));
    
    right_left_n1 = ((b(2) > c(2) && c(2) >= d(2))) || ...
        ((b(2) >= c(2) && c(2) > d(2)));
    
    
    right_left_n2 = ((a(2) > b(2) && b(2) >= c(2))) || ...
        ((a(2) >= b(2) && b(2) > c(2)));
    
    % curve from top to bottom
    % z-coordinate getting bigger
    bottom_top_h1 = ((a(3) < b(3) && b(3) <= c(3) && c(3) <= d(3))) || ...
        ((a(3) <= b(3) && b(3) < c(3) && c(3) <= d(3))) || ...
        ((a(3) <= b(3) && b(3) <= c(3) && c(3) < d(3)));
    
    bottom_top_n1 = ((b(3) < c(3) && c(3) <= d(3))) || ...
        ((b(3) <= c(3) && c(3) < d(3)));
    
    bottom_top_n2 = ((a(3) < b(3) && b(3) <= c(3))) || ...
        ((a(3) <= b(3) && b(3) < c(3)));
    
    % curve from top to bottom
    % z-coordinaten getting smaller
    top_bottom_h1 = ((a(3) > b(3) && b(3) >= c(3) && c(3) >= d(3))) || ...
        ((a(3) >= b(3) && b(3) > c(3) && c(3) >= d(3))) || ...
        ((a(3) >= b(3) && b(3) >= c(3) && c(3) > d(3)));
    
    top_bottom_n1 = ((b(3) > c(3) && c(3) >= d(3))) || ...
        ((b(3) >= c(3) && c(3) > d(3)));
    
    top_bottom_n2 = ((a(3) > b(3) && b(3) >= c(3))) || ...
        ((a(3) >= b(3) && b(3) > c(3)));
    
    % first the main condition
    left_right = left_right_h1;
    right_left = right_left_h1;
    % second the next condition
    if (~left_right && ~right_left)
        if (left_right_n1 || left_right_n2)
            left_right = true;
        elseif (right_left_n1 || right_left_n2)
            right_left = true;
        end
    end
    
    % first the main condition
    bottom_top = bottom_top_h1;
    top_bottom = top_bottom_h1;
    % second the next condition
    if (~bottom_top && ~top_bottom)
        if (bottom_top_n1 || bottom_top_n2)
            bottom_top = true;
        elseif (top_bottom_n1 || top_bottom_n2)
            top_bottom = true;
        end
    end

    if a(2) == b(2) && b(2) == c(2) && c(2) == d(2) || ...
            a(3) == b(3) && b(3) == c(3) && c(3) == d(3)
        bottom_top = 0;
        top_bottom = 0;
        left_right = 0;
        right_left = 0;
    end
    
     
    if abs(deltay_bc) > 0.02
        max_delta_offset_z = 0.00040;
    elseif abs(deltay_bc) <= 0.02 && abs(deltay_bc) > 0.01
         max_delta_offset_z = 0.00044;
    else
        max_delta_offset_z = 0.00047;
    end
    
    if abs(deltaz_bc) > 0.04
        max_delta_offset_y = 0.000049;
    elseif  abs(deltaz_bc) <= 0.04 &&  abs(deltaz_bc) > 0.03
        max_delta_offset_y = 0.000055;
    else
        max_delta_offset_y = 0.000062;
    end
    
        
    delta_offset_z = 0;
    delta_offset_y = 0;
    
    if left_right
        
        % bulge down
        if (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ... 
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ... 
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ... 
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ... 
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>steig_bc && steig_bc<=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc<steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0)
    
            delta_offset_z = -max_delta_offset_z;
        
        % bulge upwards
        elseif (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ... 
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ... 
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<steig_bc && steig_bc>=steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc>steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0)
            
            delta_offset_z = +max_delta_offset_z;
        end
    end

    if right_left
        
       % bulge down
        if (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<steig_bc && steig_bc>=steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc>steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ...
                (steig_ab>=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0)
    
            delta_offset_z = -max_delta_offset_z;
        
        % bulge upwards
        elseif (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>steig_bc && steig_bc<=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd>=0) || ...
                (steig_ab>=steig_bc && steig_bc<steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd>=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<=steig_cd && steig_ab<0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0)
    
            delta_offset_z = +max_delta_offset_z;
        end
    end

    if bottom_top
        
        % bulge to the left
        if (steig_ab>steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ... 
                (steig_ab>=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ... 
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) ||...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) ||...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<=steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0)
    
            delta_offset_y = -max_delta_offset_y;
        
        % bulge to the right
        elseif (steig_ab>steig_bc && steig_bc<=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc<steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...   
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...   
                (steig_ab<steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ...
                (steig_ab<=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0)
    
            delta_offset_y = +max_delta_offset_y;
        end
    end

    if top_bottom
        
        % bulge to the right
        if (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd<=0)
            
            delta_offset_y = +max_delta_offset_y;
        
        % bulge tot he left
        elseif (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab<=0 && steig_bc<=0 && steig_cd>=0) || ...
                (steig_ab<steig_bc && steig_bc<=steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab<=steig_bc && steig_bc<steig_cd && steig_ab>=0 && steig_bc>=0 && steig_cd>=0) || ...
                (steig_ab>steig_bc && steig_bc>=steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0) || ...
                (steig_ab>=steig_bc && steig_bc>steig_cd && steig_ab>=0 && steig_bc<=0 && steig_cd<=0)
            
            delta_offset_y = -max_delta_offset_y;
        end
    end
    
    % ---------------------
    % linear interpolation
    % ---------------------
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
    
    if (left_right ||right_left)
        offset_z = 0;
        damping = 1.0;
        % add a curvature
        for i = 1:t(end)
            
            if i <= t(end)/2 + 1
                
                damping = damping - 0.06;
                offset_z = offset_z + delta_offset_z * damping;
                 
            else
                 
                damping = damping + 0.06;
                offset_z = offset_z - delta_offset_z * damping;
                 
            end
            newpoint(3,i) = newpoint(3,i) + offset_z;
        end   
    end
    
    if (top_bottom ||bottom_top)
        offset_y = 0;
        damping = 1.0;
        % add a curvature
        for i = 1:t(end)
            
            if i <= t(end)/2 + 1
                 
                damping = damping - 0.06;
                offset_y = offset_y + delta_offset_y  * damping;
                
            else
                 
                damping = damping + 0.06;
                offset_y = offset_y - delta_offset_y  * damping;
                 
            end
            newpoint(2,i) = newpoint(2,i) + offset_y;
        end 
    end
end
