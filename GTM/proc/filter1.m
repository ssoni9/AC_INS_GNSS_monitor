%{
Filename: filter.m
by, Shivam Soni - 12/04/2021
INS monitor algorithm implementation
%}

clearvars; clc; close all;
set(0,'defaultTextInterpreter','latex');set(0,'defaultAxesFontSize',20);
set(0,'defaultLegendInterpreter','latex');
%% Load the necessary data files: 
load("sim_data/GNSS_Data.mat"); GPS_sat = GPS_sat';
t_vals = table2array(readtable("sim_data/tout.csv"));
eoms = table2struct(readtable("sim_data/AC_eom.csv"));
xdot = table2struct(readtable("sim_data/AC_xdot.csv"));
% aux_params = table2struct(readtable("sim_data/AC_aux.csv"));

%% Assign SVID's:
t_len = length(t_vals);
satIDs = 1:27;
GPS_sat(1).svid = [];
for ind = 1:t_len
    GPS_sat(ind).svid = satIDs(GPS_sat(ind).vis);
end

%% Create a struct for the problem:
prob = struct;
ft2met = 0.3048;
r2d = 180/pi;
for ind = 1:t_len
    prob(ind).t = t_vals(ind);
    prob(ind).rho = GPS_sat(:).rho;
    prob(ind).a = [xdot(ind).ub_dot; xdot(ind).vb_dot; xdot(ind).wb_dot];
    prob(ind).w = [eoms(ind).pb; eoms(ind).qb; eoms(ind).rb];
    prob(ind).sat_pos = GPS_sat(ind).pos_vis';
    prob(ind).svid = GPS_sat(ind).svid;
    prob(ind).euls = [eoms(ind).phi; eoms(ind).theta; eoms(ind).psi];
    prob(ind).LLA_GT = [eoms(ind).latitude*r2d; eoms(ind).longitude*r2d; eoms(ind).altitude*ft2met];
    prob(ind).x0_GT = lla2ecef(prob(ind).LLA_GT', 'WGS84')';
end

%% Clear structs to free-up memory: 
clear GPS_sat t_vals eoms xdot

%% Run EKF: 
% R and Q parameters:
r = 8;
q = 3;
% Initial State:
mu = [prob(1).x0_GT; 0; 0; 0; 0; 0; 0; 0];
clear prob_tight
tic
prob_tight = tightly_coupled(prob, q, r, mu);
t = toc;
fprintf('Time per estimation %f\n', t/length(prob_tight))

%% Functions: 


function prob_out = tightly_coupled(prob, q, r, mu)
    % Output Problem:
    prob_out = struct();
    
    % Initialize Covariance: 
    P = diag(rand(10,1));
    
    % Process Noise:
    Q = q*diag(ones(10,1));
    
    % State:
    prob_out(1).x0_kf = mu(1:3);
    prob_out(1).bu_kf = mu(4);
    prob_out(1).v0_kf = mu(5:6);
    prob_out(1).ang_kf = mu(7:9);
    prob_out(1).res = [];
    prob_out(1).P = [];
    prob_out(1).S0 = [];
    prob_out(1).H = [];
    prob_out(1).h = [];

    
    T_tot = length(prob);
    for ind_t = 1:T_tot       
        if ind_t == T_tot
            dt = prob(ind_t).t - prob(ind_t-1).t;
        else
            dt = prob(ind_t+1).t - prob(ind_t).t;
        end

        % Control Matrices:
        A = diag(ones(10,1)) + diag([dt*ones(3,1); 0; 0; 0], 4);
        B = [zeros(4,6); diag(dt*ones(6,1))];
        
        % Jacobian:
        meas_len = length(prob(ind_t).rho);
        [H, h] = meas_jacs(prob, prob_out, ind_t, meas_len, 1);

        % Measurement noise
        R = r*diag(ones(meas_len,1));

        % Measurement:
        z = [prob(ind_t).rho]';

        % Disturbances or Acceleration Inputs:
        acc = prob(ind_t).a;
        gyr = prob(ind_t).w;
        euls = prob(ind_t).euls;
        u_bod = [acc; gyr];
        if ind_t == 1
            u_ecef = body2ecef(u_bod, euls, prob_out(ind_t).x0_kf);
        else
            u_ecef = body2ecef(u_bod, euls, prob_out(ind_t-1).x0_kf);
        end
        
        % Call Kalman Filter
        [mu, P] = kalman_filter_ekf(A, B, H, R, Q, P, mu, z, u_ecef, h);

        % Save updated states:
        prob_out(ind_t+1).x0_kf = mu(1:3);
        prob_out(ind_t+1).bu_kf = mu(4);
        prob_out(ind_t+1).v0_kf = mu(5:7);
        prob_out(ind_t+1).ang_kf = mu(8:10);

        % Compute Residues:
        [H, h] = meas_jacs(prob, prob_out, ind_t, meas_len, 2);
        H = H(1:end, 1:4);
        S0 = (H'*H)\H'; 
        prob_out(ind_t).res = z - h;
        prob_out(ind_t).P = P;
        prob_out(ind_t).S0 = S0;
        prob_out(ind_t).H = H;
        prob_out(ind_t).h = h;
    end
end


function [mu_t_t, P_t_t] = kalman_filter_ekf(A, B, H, R, Q, P_tm_tm, mu_tm_tm, z_t, u_t, h_mu_t_tm)
    % Predict:
    mu_t_tm = A*mu_tm_tm + B*u_t;
    P_t_tm = A*P_tm_tm*A' + Q;
    
    % Update:
    y_t = z_t - h_mu_t_tm;
    K_t = P_t_tm*H'/(R + H*P_t_tm*H');
    mu_t_t = mu_t_tm + K_t*y_t;
    temp = K_t*H; temp = eye(size(temp)) - temp;
    P_t_t = temp*P_t_tm*temp' + K_t*R*K_t';
end


function u_ecef = body2ecef(u_bod, euls, pos)
    acc = reshape(u_bod(1:3),3,[]);
    gyr = reshape(u_bod(4:6),3,[]);
    pos = reshape(pos, 1, 3);
    lla = ecef2lla(pos, 'WGS84');
    R_ecef2ned = RotEcef2Ned(lla(1), lla(2));
    euls = reshape(euls, 1, []); %(end:-1:1)
    R_ned2bod = eul2rotm(euls,'ZYX');
    avec_ecef = R_ecef2ned'*R_ned2bod'*acc;
    wvec_ecef = R_ecef2ned'*R_ned2bod'*gyr;
    u_ecef = reshape([reshape(avec_ecef, 1, 3), reshape(wvec_ecef, 1, 3)],[],1);
end

function Re2n = RotEcef2Ned(latDeg, lonDeg)
    D2R = pi/180; %degrees to radians scale factor
    latRad=D2R*latDeg(:); lonRad=D2R*lonDeg(:);
    
    clat = cos(latRad);
    slat = sin(latRad);
    clon = cos(lonRad);
    slon = sin(lonRad);
    
    Re2n = zeros(3,3);
    Re2n(1,1) = -slat.*clon;
    Re2n(1,2) = -slat.*slon;
    Re2n(1,3) = clat;
    
    Re2n(2,1) = -slon;
    Re2n(2,2) = clon;
    Re2n(2,3) = 0;
    
    Re2n(3,1) = -clat.*clon;
    Re2n(3,2) = -clat.*slon;
    Re2n(3,3) = -slat;
end

function [H, h] = meas_jacs(prob, prob_out, ind_t, meas_len, flag)
    H = zeros(meas_len, 10); h = zeros(meas_len, 1);
    if flag == 1
        for ind_sat = 1:meas_len
            if ind_t == 1
                X = prob(ind_t).sat_pos(1,ind_sat) - prob_out(ind_t).x0_kf(1);
                Y = prob(ind_t).sat_pos(2,ind_sat) - prob_out(ind_t).x0_kf(2);
                Z = prob(ind_t).sat_pos(3,ind_sat) - prob_out(ind_t).x0_kf(3);
            else
                X = prob(ind_t).sat_pos(1,ind_sat) - prob_out(ind_t-1).x0_kf(1);
                Y = prob(ind_t).sat_pos(2,ind_sat) - prob_out(ind_t-1).x0_kf(2);
                Z = prob(ind_t).sat_pos(3,ind_sat) - prob_out(ind_t-1).x0_kf(3);
            end    
            eta = sqrt(X^2 + Y^2 + Z^2);
            H(ind_sat,:) = [-X/eta, -Y/eta, -Z/eta, 1, 0, 0, 0, 0, 0, 0];
            if ind_t == 1
                h(ind_sat,1) = eta + prob_out(ind_t).bu_kf;
            else
                h(ind_sat,1) = eta + prob_out(ind_t-1).bu_kf;
            end
        end
    elseif flag == 2
        for ind_sat = 1:meas_len
            X = prob(ind_t).sat_pos(1,ind_sat) - prob_out(ind_t).x0_kf(1);
            Y = prob(ind_t).sat_pos(2,ind_sat) - prob_out(ind_t).x0_kf(2);
            Z = prob(ind_t).sat_pos(3,ind_sat) - prob_out(ind_t).x0_kf(3);
            eta = sqrt(X^2 + Y^2 + Z^2);
            H(ind_sat,:) = [-X/eta, -Y/eta, -Z/eta, 1, 0, 0, 0, 0, 0, 0];
            h(ind_sat,1) = eta + prob_out(ind_t).bu_kf;
        end

    end

%     for ind_sat = 1:meas_len
%         if flag == 1
%             X = prob(ind_t).sat_pos(1,ind_sat) - prob_out(ind_t).x0_kf(1);
%             Y = prob(ind_t).sat_pos(2,ind_sat) - prob_out(ind_t).x0_kf(2);
%             Z = prob(ind_t).sat_pos(3,ind_sat) - prob_out(ind_t).x0_kf(3);   
%             eta = sqrt(X^2 + Y^2 + Z^2);
%             H(ind_sat,:) = [-X/eta, -Y/eta, -Z/eta, 1, 0, 0, 0, 0, 0, 0];
%             if ind_t == 1
%                 h(ind_sat,1) = eta + prob_out(ind_t).bu_kf;
%             else
%                 h(ind_sat,1) = eta + prob_out(ind_t-1).bu_kf;
%             end
%         elseif flag == 2
%              X = prob(ind_t).sat_pos(1,ind_sat) - prob_out(ind_t).x0_kf(1);
%             Y = prob(ind_t).sat_pos(2,ind_sat) - prob_out(ind_t).x0_kf(2);
%             Z = prob(ind_t).sat_pos(3,ind_sat) - prob_out(ind_t).x0_kf(3);   
%             eta = sqrt(X^2 + Y^2 + Z^2);
%             H(ind_sat,:) = [-X/eta, -Y/eta, -Z/eta, 1, 0, 0, 0, 0, 0, 0];
%             if ind_t == 1
%                 h(ind_sat,1) = eta + prob_out(ind_t).bu_kf;
%             else
%                 h(ind_sat,1) = eta + prob_out(ind_t-1).bu_kf;
%             end
%         end
%     end
end