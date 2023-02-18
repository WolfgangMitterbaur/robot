% Trajectory generation for polynomial profiles
% Wolfgang Mitterbaur

%% initialization of the robot
clear, clc, close all
load gen3              % load the robot
load exampleHelperKINOVAGen3GripperCollRRT.mat
load gen3positions     % load the initial position
endEffector = 'Gripper';    % definition of the end-effector
numJoints = numel(gen3.homeConfiguration);

% load the waypoints
CalculateWaypoints;
waypointTimes = 0:0.5:47.5; % all points %  0:0.5:9; % only the "h"
% trajectory sample time
ts = 0.05;
trajTimes = 0:ts:waypointTimes(end);
% Acceleration (quintic only)
waypointAccels = zeros(size(waypointVels));

% Acceleration times (trapezoidal only)
waypointAccelTimes = diff(waypointTimes)/2;


% frame and portal
ax = gca;
plane = collisionBox(1.0, 1.4, 0.02);
plane.Pose = trvec2tform([0.0 0.0 0.02]);
boardtop = collisionBox(0.02, 1.4, 0.02);
boardtop.Pose = trvec2tform([0.55 0 0.9]);
boardbottom = collisionBox(0.02, 1.4, 0.02);
boardbottom.Pose = trvec2tform([0.55 0 0.18]);
boardleft = collisionBox(0.02, 0.02, 0.74);
boardleft.Pose = trvec2tform([0.55 -0.7 0.54]);
boardright = collisionBox(0.02, 0.02, 0.74);
boardright.Pose = trvec2tform([0.55 0.7 0.54]);

endEffectorPose = getTransform(gen3, jointAnglesHome', "Gripper");
pencil = collisionBox(0.15, 0.02, 0.02);   % the pencil for the robot
pencil_pos = endEffectorPose(1:3,4)';
pencil.Pose = trvec2tform(pencil_pos);

ik = inverseKinematics('RigidBodyTree', gen3);
ikRobotWeights = [1 1 1 1 1 1]; % definition of weights of each joint
ikInitConfig = gen3.homeConfiguration;

% --- show the robot in home configuration
show(robot, jointAnglesHome', 'Frames', 'off', 'PreservePlot', false);
xlim([-1 1]), ylim([-1 1]), zlim([0 1.2])

hold on
% --- plot the waypints (trajectory)
% hTraj = plot3(waypoints(1,1), waypoints(2,1), waypoints(3,1),'b.-');
% --- plot a linear connection line between the waypoints
% plot3(waypoints(1,:), waypoints(2,:), waypoints(3,:), 'ro', 'LineWidth', 0.001);
% --- plot a linear connection line between the waypoints
% set(hTraj, 'xdata', waypoints(1,:), 'ydata', waypoints(2,:), 'zdata', waypoints(3,:));


%% Generate trajectory
% --- trapezodial velocity profiles
% [q, qd, qdd] = trapveltraj(waypoints, numel(trajTimes), ...
%             'AccelTime',repmat(waypointAccelTimes,[3 1]), ... 
%             'EndTime',repmat(diff(waypointTimes),[3 1]));
% --- cubic polynomial
% [q, qd, qdd] = cubicpolytraj(waypoints, waypointTimes, trajTimes, ... 
%            'VelocityBoundaryCondition', waypointVels);
% --- quintic polynomial
% [q, qd, qdd] = quinticpolytraj(waypoints, waypointTimes, trajTimes, ... 
%             'VelocityBoundaryCondition', waypointVels, ...
%             'AccelerationBoundaryCondition', waypointAccels);
% --- spline polynomial
[q, qd, qdd] = bsplinepolytraj(waypoints, waypointTimes([1 end]), trajTimes);

% --- plot all calculated points
%set(hTraj, 'xdata', q(1,:), 'ydata', q(2,:), 'zdata', q(3,:));
% --- plot all calculated points with a cooridnate frame
% plotTransforms(q', repmat([1 0 0 0], [size(q,2) 1]), 'FrameSize', 0.05);
% --- plot the trajectories as a function of time
% PlotTrajectory(trajTimes, q, qd, qdd, 'Names', ["X","Y","Z"], 'WaypointTimes', waypointTimes)

% plot the frame and the portal
[~, patchObj] = show(plane,'Parent', ax);
patchObj.FaceColor = [0.53 0.53 0.53];
[~, patchObj] = show(boardtop,'Parent', ax);
patchObj.FaceColor = [0.79 0.88 1.0];
[~, patchObj] = show(boardbottom,'Parent', ax);
patchObj.FaceColor = [0.79 0.88 1.0];
[~, patchObj] = show(boardleft,'Parent', ax);
patchObj.FaceColor = [0.79 0.88 1.0];
[~, patchObj] = show(boardright,'Parent', ax);
patchObj.FaceColor = [0.79 0.88 1.0];

% inverse kinematic for all waypoints
for index = 1:numel(trajTimes) 
    % calculate inverse kinematics
    robotPose = trvec2tform(q(:,index)');
    [config, info] = ik(endEffector, robotPose, ikRobotWeights, ikInitConfig);

    % --- show the robot
    show(robot, config, 'Frames', 'off', 'PreservePlot', false);
    title(['Trajectory at t = ' num2str(trajTimes(index))])
    % --- plot the trajectory points
    plot3(q(1, index), q(2, index), q(3, index), 'b.-');

    drawnow    
end

pencil = collisionBox(0.1, 0.02, 0.02);   % the pencil for the robot
endEffectorPose = getTransform(robot, config, "gripper");
pencil_pos = endEffectorPose(1:3,4)';
pencil.Pose = trvec2tform(pencil_pos);
show(pencil);
