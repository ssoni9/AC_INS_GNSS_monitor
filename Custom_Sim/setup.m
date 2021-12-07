clearvars; clc; close all; 
open("ACsim.slx")
% Link: https://www.mathworks.com/help/satcom/ug/estimate-gnss-receiver-position.html


%%
 Fs = 1;
numSamples = 1000;
t = 0:1/Fs:(numSamples-1)/Fs;
% LLA position for Natick, MA
refLocNatick = [42.2825 -71.343 53.0352];

gps = gpsSensor('SampleRate', Fs, ...
    'ReferenceLocation', refLocNatick);

pos = zeros(numSamples, 3);
vel = zeros(numSamples, 3);

llaMeas = gps(pos, vel);

subplot(3, 1, 1)
plot(t, llaMeas(:,1))
title('Latitude')
xlabel('s')
ylabel('degrees')

subplot(3, 1, 2)
plot(t, llaMeas(:,2))
title('Longitude')
xlabel('s')
ylabel('degrees')

subplot(3, 1, 3)
plot(t, llaMeas(:,3))
title('Altitude')
xlabel('s')
ylabel('m')