% Trajectory Generation
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
waypointTimes = 0:0.5:47.5; % all points 
%waypointTimes =  0:0.5:9; % only the "h"
% trajectory sample time
ts = 0.04546;       % 47.5 / 1045 (totalTime / number of points)
%ts = 0.023810;      % 47.5 / 1995 (totalTime / number of points)
%ts = 0.02525; %0.048129;      % 9 / 187(357) (totalTime / number of points)
trajTimes = 0:ts:waypointTimes(end);
% acceleration times 
% waypointAccelTimes = diff(waypointTimes)/2;

% process the interpolation
ProcessInterpos;

% frame and portal
ax = gca;
plane = collisionBox(1.0, 1.4, 0.02);   % the portal for the robot
plane.Pose = trvec2tform([0.0 0.0 0.02]);
boardtop = collisionBox(0.02, 1.4, 0.02);   % the top part of the frame
boardtop.Pose = trvec2tform([0.63 0 0.9]);
boardbottom = collisionBox(0.02, 1.4, 0.02);    % the bottom part of the frame
boardbottom.Pose = trvec2tform([0.63 0 0.18]);
boardleft = collisionBox(0.02, 0.02, 0.74);     % the left part of the frame
boardleft.Pose = trvec2tform([0.63 -0.7 0.54]);
boardright = collisionBox(0.02, 0.02, 0.74);    % the right part of the frame
boardright.Pose = trvec2tform([0.63 0.7 0.54]);

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
%plotpoints = waypoints;
%plotpoints = plotpoints + [0.07 0 0]';
%hTraj = plot3(plotpoints(1,1), plotpoints(2,1), plotpoints(3,1),'b.-');
% --- plot a linear connection line between the waypoints
%plot3(plotpoints(1,:), plotpoints(2,:), plotpoints(3,:), 'ro', 'LineWidth', 0.001);
% --- plot a linear connection line between the waypoints
% set(hTraj, 'xdata', waypoints(1,:), 'ydata', waypoints(2,:), 'zdata', waypoints(3,:));


%% Generate trajectory
% --- plot all calculated points
% set(hTraj, 'xdata', q(1,:), 'ydata', q(2,:), 'zdata', q(3,:));
% --- plot all calculated points with a cooridnate frame
% plotTransforms(q', repmat([1 0 0 0], [size(q,2) 1]), 'FrameSize', 0.05);
% --- plot the trajectories as a function of time
%PlotTrajectory2(trajTimes, q, qd, qdd, qd_f, qdd_f, 'Names', ["X","Y","Z"], 'WaypointTimes', waypointTimes)

%plot the frame and the portal
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
% [~, patchObj] = show(pencil,'Parent', ax);
% patchObj.FaceColor = [1.0 0.0 0.0];

plotq = q  + [0.07 0 0]';
qs = zeros(numel(trajTimes), 7);

% inverse kinematic for all waypoints
for index = 1:numel(trajTimes) 
    % calculate inverse kinematics
    robotPose = trvec2tform(q(:,index)');
    [config, info] = ik(endEffector, robotPose, ikRobotWeights, ikInitConfig);
   
    % --- show the robot
    show(robot, config, 'Frames', 'off', 'PreservePlot', false);
  
%     hold on
% 
%     [~, patchObj] = show(plane,'Parent', ax);
%     patchObj.FaceColor = [0.53 0.53 0.53];
%     [~, patchObj] = show(boardtop,'Parent', ax);
%     patchObj.FaceColor = [0.79 0.88 1.0];
%     [~, patchObj] = show(boardbottom,'Parent', ax);
%     patchObj.FaceColor = [0.79 0.88 1.0];
%     [~, patchObj] = show(boardleft,'Parent', ax);
%     patchObj.FaceColor = [0.79 0.88 1.0];
%     [~, patchObj] = show(boardright,'Parent', ax);
%     patchObj.FaceColor = [0.79 0.88 1.0];
   
    % show the pencil
%     pencil = collisionBox(0.1, 0.02, 0.02);   % the pencil for the robot
%     endEffectorPose = getTransform(robot, config, "gripper");
%     pencil_pos = endEffectorPose(1:3,4)';
%     pencil.Pose = trvec2tform(pencil_pos);
%     show(pencil);
    
    %title(['Trajectory at t = ' num2str(trajTimes(index))])
    % --- plot the trajectory points
    plot3(plotq(1, index), plotq(2, index), plotq(3, index), 'b.-');
    
    %hold off;
    
    % store the configuration
    qs(index,:) = config(1:7); 
    drawnow    
end

pencil = collisionBox(0.1, 0.02, 0.02);   % the pencil for the robot
endEffectorPose = getTransform(robot, config, "gripper");
pencil_pos = endEffectorPose(1:3,4)';
pencil.Pose = trvec2tform(pencil_pos);
show(pencil);

% --- plot the joints
%PlotJoints(trajTimes, qs, ts, waypointTimes);

