% function to plot the position, velocity and acceleration
% Wolfgang Mitterbaur

%% function PlotTrajectory
% input values:
% t: time series
% q: position series
% qd: velocity series
% qdd: acceleration series
% qd_f: filtered velocity series
% qdd_f: filtered acceleration series
% varargin: additional input parameters

function PlotTrajectory2(t, q, qd, qdd, qd_f, qdd_f, varargin)

    order = 3; % plot position, velocity and acceleration
    
    % initialize default plot options
    numCoords = size(q,1);
    configNames = "Coordinate " + varargin{2};
      
    %waypointTimes = [];
    waypointTimes = varargin{4};

    % loop through all the coordinates and plot
    for index = 1:numCoords
        figure;
        
        % always plot time and position only
        subplot(order, 1, 1), hold on;
        plot(t, q(index,:));
        plotTimeLines(waypointTimes);
        ylabel('Position');
        title("Trajectory for " + configNames(index))
        title("" + configNames(index))
        
        subplot(order, 1, 2), hold on;
        %qd = varargin{1};
        plot(t, qd(index,:));
        t_neu = t-2/24;
        qd_f_neu = qd_f(index,:);
        plot(t_neu(6:numel(t)),qd_f_neu((6:numel(t))));
        plotTimeLines(waypointTimes);
        ylabel('Velocity');
           
        subplot(order, 1, 3), hold on;
        %qdd = varargin{2};
        plot(t, qdd(index,:));
        qdd_f_neu = qdd_f(index,:);
        plot(t_neu(6:numel(t)), qdd_f_neu(6:numel(t)));
        plotTimeLines(waypointTimes);
        ylabel('Acceleration');
 
        % label the time axis
        xlabel('Time');
       
    end

end

% function to plot vertical lines for each waypoint
function plotTimeLines(t)
    for index = 1:numel(t)
       xline(t(index), 'r--'); 
    end
end

