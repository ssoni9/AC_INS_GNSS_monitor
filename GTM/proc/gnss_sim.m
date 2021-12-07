%% Simulate GNSS measurements giving G.T. PVT of A/C:
clearvars; clc; close all;

%% Import the GT csv files:
AC_eom = readtable("sim_data/AC_eom.csv");
tvals_sim = readtable("sim_data/tout.csv");
tvals_sim = table2array(tvals_sim);

%% Get Satellite positions and velocities:

t_init = datetime("now","TimeZone","Local");
t_vec = t_init + seconds(tvals_sim);

% Preallocate:
t_len = length(tvals_sim);
GPS_sat = struct;
GPS_sat(1).pos = [];
GPS_sat(1).vel = [];
for ind = 1:t_len
    [GPS_sat(ind).pos, GPS_sat(ind).vel] = gnssconstellation(t_vec(ind));
end
%% Mask Satellties below 5 degrees of elevation:
ft2met = 0.3048;
GT_LLA_pos = [AC_eom.latitude*180/pi, AC_eom.longitude*180/pi, AC_eom.altitude*ft2met];
mask_angle = 5;
GPS_sat(1).vis = true;
GPS_sat(1).pos_vis = [];
GPS_sat(1).vel_vis = [];
for ind = 1:t_len
    rxpos = GT_LLA_pos(ind,:);
    sat_pos = GPS_sat(ind).pos;
    n_sat = size(sat_pos,1);
    vis_ind = 1;
    for sat_ind = 1:n_sat
        [~, ~, vis_temp] = lookangles(rxpos, sat_pos(sat_ind,:), mask_angle);
        GPS_sat(ind).vis(sat_ind) = vis_temp;
        if vis_temp 
            GPS_sat(ind).pos_vis(vis_ind,:) = GPS_sat(ind).pos(sat_ind,:);
            GPS_sat(ind).vel_vis(vis_ind,:) = GPS_sat(ind).vel(sat_ind,:);
            vis_ind = vis_ind + 1;
        end
    end
end

%% Get pseudoranges:
% Preallocate:
GPS_sat(1).rho = [];

for ind = 1:t_len
    rxpos = GT_LLA_pos(ind,:);
    sat_pos = GPS_sat(ind).pos_vis;
    n_sat = size(sat_pos,1);
    for sat_ind = 1:n_sat
        GPS_sat(ind).rho(sat_ind) = pseudoranges(rxpos, sat_pos(sat_ind, :), 'RangeAccuracy', 1);
    end
end

%% Save the "GPS_sat" struct:

save("sim_data/GNSS_Data2.mat",'GPS_sat')




