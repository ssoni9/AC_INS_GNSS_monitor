%% Plot LLA
% x_lla = zeros(length(prob_tight),3); x_lla_GT = x_lla;
set(0,'defaultTextInterpreter','latex'); 
set(0,'defaultAxesFontSize',20); 
cmat = parula(length(prob_tight));
for ind=1:length(prob_tight)-1
    x_lla(ind,:) = ecef2lla(prob_tight(ind).x0_kf','WGS84'); 
    x_lla_GT(ind,:) = prob(ind).LLA_GT;
end
% figure;
% for i=1:100:length(prob_tight)
%     geoplot(x_lla(i,1), x_lla(i,2),'.','Color',cmat(i,:));hold on;
% end
% title('Receiver position in LLA frame using WGS84 datum');
% colorbar('southoutside','TickLabelInterpreter','latex','FontSize',24,...
%     'TicksMode','manual','Ticks',[0, 1], 'TickLabels',{'$t = 0$', '$t = t_{end}$'})
% figure;
% for i=1:100:length(prob_tight)
%     geoplot(x_lla_GT(i,1), x_lla_GT(i,2), 'Color','r')
% end
%%
figure; hold on; 
plot(x_lla_GT(:,1), x_lla_GT(:,2),'b','LineWidth',2)
hold on; plot(x_lla(:,1), x_lla(:,2),':r'); grid minor
% hold on; plot(x_lla(:,1), x_lla(:,2),':k'); grid minor


xlabel('Latitude [deg]'); ylabel('Longitude [deg]'); legend('Ground Truth', 'EKF Estimate');
title('Multiple Faulty Satellites Trajectory')
%% 
for ind = 1:length(prob_tight)-1
    res(ind) = norm(prob_tight(ind).res);
end
hold on; plot([prob.t], 40*ones(length([prob.t]),1), 'b', 'LineWidth',2)
hold on; plot([prob.t], res,'r'); grid minor
axis([0 60 0 45])
xlabel('Time [s]'); ylabel('Residual Statistic');
legend('Error', 'Alert limit'); title('Multiple Faulty Satellites - AL Exceeded')
%%
for ind = 1:length(prob_tight)-1
    temp = (eye(7) - prob_tight(ind).h)*prob(ind).rho';
    res(ind) = sqrt(temp'*temp);
end
plot(res)

%% 
HH = prob_tight(1).H;
S0 = inv(HH'*HH)*HH';
z = prob_tight(1).h;
z = prob(1).rho';
r  = (eye(7) - HH*S0)*z
q = sqrt(r'*r)
