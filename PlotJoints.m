% Plot the joint trajectories
% Wolfgang Mitterbaur

%% calulate the positions, velocities and acceleration
% input values:
% trajTimes: the time values
% qs: the path including positions
% ts: time step
% waypointTimes: the timing law including all time steps

function PlotJoints(trajTimes, qs, ts, waypointTimes)

    % angle in grad
    qs_deg = qs*180/pi;     
    
    % filter
    h = [1/2 1/2];
    binomialCoeff = conv(h,h);
    for n = 1:4
        binomialCoeff = conv(binomialCoeff, h);
    end
    
    % velocity in grad / s
    vel = diff(qs_deg)/ts;     
    vel(1,:) = 0;
    vel(end+1,:) = 0;
    
    % filter for velocity
    vel_f = filter(binomialCoeff, 1, vel);
    
    % acceleration in grad / s2
    acc = diff(vel_f)/ts;    
    acc(1,:) = 0;
    acc(end+1,:) = 0;
    
    % filter for acceleration
    acc_f = filter(binomialCoeff, 1, acc);
    
    plotqs = qs_deg(:,1:3)';
    plotvel = vel(:,1:3)';
    plotacc = acc(:,1:3)';
    plotvel_f = vel_f(:,1:3)';
    plotacc_f = acc_f(:,1:3)';
    PlotTrajectory(trajTimes,plotqs,plotvel,plotacc,plotvel_f,plotacc_f,'Names',["Joint 1","Joint 2","Joint 3"],'WaypointTimes',waypointTimes)
    
    plotqs = qs_deg(:,4:6)';
    plotvel = vel(:,4:6)';
    plotacc = acc(:,4:6)';
    plotvel_f = vel_f(:,4:6)';
    plotacc_f = acc_f(:,4:6)';
    PlotTrajectory(trajTimes,plotqs,plotvel,plotacc,plotvel_f,plotacc_f,'Names',["Joint 4","Joint 5","Joint 6"],'WaypointTimes',waypointTimes)
    
    plotqs = qs_deg(:,7)';
    plotvel = vel(:,6:7)';
    plotacc = acc(:,6:7)';
    plotvel_f = vel_f(:,6:7)';
    plotacc_f = acc_f(:,6:7)';
    PlotTrajectory(trajTimes,plotqs,plotvel,plotacc,plotvel_f,plotacc_f,'Names',["Joint 7"],'WaypointTimes',waypointTimes)

end


