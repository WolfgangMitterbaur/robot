% Function to flatten the curve
% Wolfgang Mitterbaur

%% Flatten the curve

no = numel(q(2,:))-3;
delta_y_12 = zeros(1, no); 
delta_y_23 = zeros(1, no); 
delta_y_34 = zeros(1, no); 
delta_z_12 = zeros(1, no); 
delta_z_23 = zeros(1, no); 
delta_z_34 = zeros(1, no); 

for i = 1:no 
    delta_y_12(i) = q(2, i+1) - q(2, i);
    delta_y_23(i) = q(2, i+2) - q(2, i+1);
    delta_y_34(i) = q(2, i+3) - q(2, i+2);

    if delta_y_23(i) > 0 && delta_y_12(i) > 0 && delta_y_34(i) > 0 && ...
        q(2, i) < q(2, i+1) && q(2, i+1) < q(2, i+2) && q(2, i+2) < q(2, i+3) && ...
        delta_y_23(i) > delta_y_12(i)*1.15 && delta_y_23(i) > delta_y_34(i)*1.15 && ...
        (( q(2, i) < 0 && q(2, i+1) < 0 && q(2, i+2) < 0 && q(2, i+3) < 0) || ...
         (q(2, i) > 0 && q(2, i+1) > 0 && q(2, i+2) > 0 && q(2, i+3) > 0))
            q(2, i+2) = q(2, i+2) - delta_y_23(i)* 0.25;
    end

    if delta_y_23(i) < 0 && delta_y_12(i) < 0 && delta_y_34(i) < 0  && ...
        q(2, i) > q(2, i+1) && q(2, i+1) > q(2, i+2) && q(2, i+2) > q(2, i+3) && ...
        delta_y_23(i) < delta_y_12(i)*1.15 && delta_y_23(i) < delta_y_34(i)*1.15 && ...
        (( q(2, i) < 0 && q(2, i+1) < 0 && q(2, i+2) < 0 && q(2, i+3) < 0) || ...
         (q(2, i) > 0 && q(2, i+1) > 0 && q(2, i+2) > 0 && q(2, i+3) > 0))
            q(2, i+2) = q(2, i+2) - delta_y_23(i)* 0.25; 
    end

    delta_z_12(i) = q(3, i+1) - q(3, i);
    delta_z_23(i) = q(3, i+2) - q(3, i+1);
    delta_z_34(i) = q(3, i+3) - q(3, i+2);

    if delta_z_23(i) > 0 && delta_z_12(i) > 0 && delta_z_34(i) > 0 && ...
        q(3, i) < q(3, i+1) && q(3, i+1) < q(3, i+2) && q(3, i+2) < q(3, i+3) && ...
        delta_z_23(i) > delta_z_12(i)*1.15 && delta_z_23(i) > delta_z_34(i)*1.15 && ...
        (( q(3, i) < 0 && q(3, i+1) < 0 && q(3, i+2) < 0 && q(3, i+3) < 0) || ...
         (q(3, i) > 0 && q(3, i+1) > 0 && q(3, i+2) > 0 && q(3, i+3) > 0))
            q(3, i+2) = q(3, i+2) - delta_z_23(i)* 0.25;
    end
     
    if delta_z_23(i) < 0 && delta_z_12(i) < 0 && delta_z_34(i) < 0  && ...
        q(3, i) > q(3, i+1) && q(3, i+1) > q(3, i+2) && q(3, i+2) > q(3, i+3) && ...
        delta_z_23(i) < delta_z_12(i)*1.15 && delta_z_23(i) < delta_z_34(i)*1.15 && ...
        ((q(3, i) < 0 && q(3, i+1) < 0 && q(3, i+2) < 0 && q(3, i+3) < 0) || ...
         (q(3, i) > 0 && q(3, i+1) > 0 && q(3, i+2) > 0 && q(3, i+3) > 0))
            q(3, i+2) = q(3, i+2) - delta_z_23(i)* 0.25;
    end

end









gradient = zeros(1, no); 
delta_gradient = zeros(1, no); 
for i = 1:no 

    delta_y = q(2, i+1) - q(2, i);
    delta_z = q(3, i+1) - q(3, i);
    
%     if delta_y == 0 
%         gradient(i) = 100;
%     else
        gradient(i) = delta_z / delta_y;
%     end

    if gradient(i) < -50
        gradient(i) = -50;
    end
     if gradient(i) > 50
        gradient(i) = 50;
    end

   

end

% figure
% plot(gradient)
% %hold on

for i = 1:no-1 
 delta_gradient(i) =  gradient(i+1) - gradient(i);


    if delta_gradient(i) < -50
        delta_gradient(i) = -50;
    end
    if delta_gradient(i) > 50
        delta_gradient(i) = 50;
    end

end

% figure
% plot(delta_gradient)

% trajTimes1 = trajTimes(1:354);
% PlotTrajectory(trajTimes1, delta_gradient, qd, qdd, qd_f, qdd_f, 'Names', ["delta gradient","Y","Z"], 'WaypointTimes', waypointTimes)
% PlotTrajectory(trajTimes1, gradient, qd, qdd, qd_f, qdd_f, 'Names', ["gradient","Y","Z"], 'WaypointTimes', waypointTimes)


