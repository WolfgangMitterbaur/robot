%velo_unfiltered = qd(2,:);
% coeffvelo = ones(1, 10)/10;
% velo_filtered = filter(coeffvelo, 1, velo_unfiltered);
% fDelay_velo = (length(coeffvelo)-1)/2;
% 
% figure;
% plot (trajTimes, velo_unfiltered);
% hold on;
% plot (trajTimes - fDelay_velo/10, velo_filtered);


% accel_unfiltered = qd(3,:);
% coeffaccel = ones(1, 10)/10;
% accel_filtered = filter(coeffaccel, 1, accel_unfiltered);
% fDelay_accel = (length(coeffaccel)-1)/2;
% 
% figure;
% plot (trajTimes, accel_unfiltered);
% hold on;
% plot (trajTimes- fDelay_accel/10, accel_filtered);


figure;
velo_unfiltered = qd(2,:);
h = [1/2 1/2];
binomialCoeff = conv(h,h);
for n = 1:4
    binomialCoeff = conv(binomialCoeff,h);
end
fDelay = (length(binomialCoeff)-1)/2;
binomialMA = filter(binomialCoeff, 1, velo_unfiltered);
plot(trajTimes, velo_unfiltered, trajTimes - fDelay / 24, binomialMA)

% figure
% alpha = 0.45;
% exponentialMA = filter(alpha, [1 alpha-1], velo_unfiltered);
% plot(trajTimes, velo_unfiltered, ...
%      trajTimes-fDelay/24,binomialMA, ...
%      trajTimes-1/24,exponentialMA)