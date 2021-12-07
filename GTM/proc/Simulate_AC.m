%% Simulate:

% Load nominal starting parameter set
MWS_Nominal=init_design();

% Set nominal trim to climbing turn
MWS_Nominal.DamageCase=0;
MWS_Nominal.DamageOnsetTime=120;
loadmws(MWS_Nominal);
[MWS_Nominal,Xtrim,Fcond,Err]=trimgtm(struct('eas',100, 'pitch',0,...
                                             'yawrate',0,'roll',0,...
                                             'rollrate',0, 'pitchrate',0,...
                                             'gndtrack', 180));
% Fcond.pitch = -3;
% Fcond.altitude = 1500;

% Simulate flight without damage
loadmws(MWS_Nominal);
fprintf(1,'Simulating...');
sim('gtm_design',[0 60]);
fprintf(1,'Done\n');
tout1=tout;Lon1=sout.eom.longitude*180/pi; Lat1=sout.eom.latitude*180/pi; Alt1=sout.eom.altitude; 
eas1=sout.aux.eas;alpha1=sout.aux.alpha; beta1=sout.aux.beta; gamma1=sout.aux.gamma; 
pqr1=180/pi*[sout.eom.pb sout.eom.qb sout.eom.rb];
phi1=180/pi*sout.eom.phi;
theta1=180/pi*sout.eom.theta;
psi1=180/pi*sout.eom.psi;

%% Plots:
set(0,'defaultTextInterpreter','latex'); 
set(0,'defaultAxesFontSize',20); 

h1=plot3(Lat1,Lon1,Alt1','b'); grid on,hold on
xlabel('Latitude (deg)'),ylabel('Longitude (deg)'),zlabel('Altitude (ft)');

%% Save the required simulation outputs into ".csv" files:
% Save time: 
t_tab = array2table(tout, 'VariableNames', {'t'});
writetable(t_tab, 'sim_data/tout.csv')
% Save state:
state_table = struct2table(sout.eom);
writetable(state_table, 'sim_data/AC_eom.csv')

% Save state derivative (acceleration)
state_table = struct2table(sout.xdot);
writetable(state_table, 'sim_data/AC_xdot.csv')

% Save AC Params
state_table = struct2table(sout.AC_Params);
writetable(state_table, 'sim_data/AC_Params.csv')

% Save aux Params
struct_new = sout.aux;
struct_new = rmfield(struct_new, 'DCM');
state_table = struct2table(struct_new);
writetable(state_table, 'sim_data/AC_aux.csv')
writematrix(sout.aux.DCM, 'sim_data/AC_DCM.csv')
% Read DCM as: 
% XX = readmatrix('sim_data/AC_DCM.csv');
% XX = reshape(XX, 3, 3,[]);

%% Extra Stuff:
%{
%% Save LLA as .csv:
LLA_tab = array2table([tout1, Lat1, Lon1, Alt1]);
LLA_tab.Properties.VariableNames = [{'t'}, {'Lat'}, {'Lon'}, {'Alt'}];
writetable(LLA_tab, 'processing/AC_LLA.csv')

%%
baselineApproachTrajectory = geoTrajectory([37.2830 -121.8395 2000; ...
    37.3630 -121.9287 62],[0;400], ...
    'ClimbRate', [6; 3.75]);
viewer = helperPertScenarioGlobeViewer;
positionCamera(viewer, [42.3072 -70.8463 12455], [0 -34 335]);
plotTrajectory(viewer, baselineApproachTrajectory, 'Color', [15 255 255]/255, "Width", 1);
%}
