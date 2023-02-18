% Process the interpolation between all waypoints
% Wolfgang Mitterbaur

%% execute all waypoints
no_points = 97 - 4;     % all points
%no_points = 19 - 4;     % only point for letter h

% number of interpolation points between 2 waypoints
no_interpol = 10;   
% total number of points
no_way = 95 * no_interpol + 95;      %% all points
%no_way = 17 * no_interpol + 17;      %% all points

% array with new waypoints
q = zeros(3, no_way);  
% array with velocities
qd = zeros(3, no_way); 
% array with acceleration
qdd = zeros(3, no_way);

% the index of the inital array of waypoints
initial_pos = 1;
% the index of the final array of waypoints
curr_pos = 1;       

%% the first point
a = waypoints(:,initial_pos);
b = waypoints(:,initial_pos + 1);
c = waypoints(:,initial_pos + 2);
q(:, curr_pos) = a;     % add the existing waypoint
curr_pos = curr_pos + 1;

newpoint = Interpol(a, a, b, c, no_interpol);
for i = 1:no_interpol 
    q(:, i + curr_pos - 1) = newpoint(:, i);    % add all new waypoints
end
curr_pos = curr_pos + no_interpol;

%% all points between the first and the last point
for j = 1:no_points 
    a = waypoints(:,initial_pos);
    b = waypoints(:,initial_pos + 1);
    c = waypoints(:,initial_pos + 2);
    d = waypoints(:,initial_pos + 3);
    
    q(:, curr_pos) = b;     % add the existing waypoint
    curr_pos = curr_pos + 1;
    
    newpoint = Interpol(a, b, c, d, no_interpol);

    for i = 1:no_interpol 
        q(:, i + curr_pos-1) = newpoint(:, i);     % add all new waypoints
    end
    curr_pos = curr_pos + no_interpol;
    initial_pos = initial_pos + 1; % next character

    % 5 points at the connection waypoint
    change_point = curr_pos - no_interpol - 1;   
    first_point = change_point - 3;
    last_point = change_point + 3;
    a = q(:, first_point);
    b = q(:, last_point);
    c = q(:, change_point);
    
    % check if the points are rising up or rising down
    if a(2) < c(2) &&  c(2) < b(2) || a(3) < c(3) && c(3) < b(3) || ...
        a(2) > c(2) &&  c(2) > b(2) || a(3) > c(3) && c(3) > b(3)
        newpoint = Interpol2(a, b, c, 5);
        for i = 1:5
            q(:, first_point + i) = newpoint(:, i);
        end
    end
end

%% last point
a = waypoints(:, initial_pos);
b = waypoints(:, initial_pos + 1);
c = waypoints(:, initial_pos + 2);
q(:, curr_pos) = b;         % add the existing waypoint
curr_pos = curr_pos + 1;

newpoint = Interpol(a, b, c, c, 20);
for i = 1:no_interpol 
    q(:, i + curr_pos-1) = newpoint(:, i);  % add all new waypoints
end
curr_pos = curr_pos + no_interpol;

% 5 points at the connection point
change_point = curr_pos - no_interpol-1;   
last_point = change_point + 3;
first_point = change_point - 3;
a = q(:, first_point);
b = q(:, last_point);
c = q(:, change_point);

newpoint = Interpol2(a, b, c, 5);
% check if the points are rising up or rising down
if a(2) < c(2) &&  c(2) < b(2) || a(3) < c(3) && c(3) < b(3) || ...
    a(2) > c(2) &&  c(2) > b(2) || a(3) > c(3) && c(3) > b(3)
    newpoint = Interpol2(a, b, c, 5);
    for i = 1:5
        q(:, first_point + i) = newpoint(:, i);
    end
end


%% calculate velocity and acceleration with differences

qd(:,1) = [0, 0, 0];

for i = 2:no_way-1 
    qd(:,i) = (q(:,i) - q(:,i-1)) / ts;
 end

% filter paramameter
h = [1/2 1/2];
binomialCoeff = conv(h,h);
for n = 1:4
    binomialCoeff = conv(binomialCoeff, h);
end
fDelay = (length(binomialCoeff)-1)/2;

% filter for velocity
velo_unfiltered = qd(2,:);
v_binomialMA = filter(binomialCoeff, 1, velo_unfiltered);
% figure;
% plot(trajTimes, velo_unfiltered, trajTimes - fDelay / 24, v_binomialMA)
% title("velocity 2");
qd_f(2,:) = v_binomialMA;

velo_unfiltered = qd(3,:);
v_binomialMA = filter(binomialCoeff, 1, velo_unfiltered);
% figure;
% plot(trajTimes, velo_unfiltered, trajTimes - fDelay / 24, v_binomialMA)
% title("velocity 3");
qd_f(3,:) = v_binomialMA;

qdd(:,1) = [0, 0, 0];
for i = 2:no_way-1 
    qdd(:,i) = (qd_f(:,i)- qd_f(:,i-1)) / ts;
end
% filter for acceleration
accel_unfiltered = qdd(2,:);
a_binomialMA = filter(binomialCoeff, 1, accel_unfiltered);
% figure;
% plot(trajTimes, accel_unfiltered, trajTimes - fDelay / 24, a_binomialMA)
% title("acceleration 2");
qdd_f(2,:) = a_binomialMA;

accel_unfiltered = qdd(3,:);
a_binomialMA = filter(binomialCoeff, 1, accel_unfiltered);
% figure;
% plot(trajTimes, accel_unfiltered, trajTimes - fDelay / 24, a_binomialMA)
% title("acceleration 3");
qdd_f(3,:) = a_binomialMA;

%% Flatten the curve
% flatten the curve of interpolation points
Flatten;
