function plot_six(U_true_c1,U_sim_c1,U_est_c1,U_true_c2,U_sim_c2,U_est_c2,model_true)
res = 600;
fontsize = 20;
% true
figure;
subplot(231)
plot_sc(U_true_c1,model_true);
% title('True density','interpreter', 'latex');
title('True denstiy $\rho_1$','interpreter', 'latex');
set(gca,'fontsize',fontsize)
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[100  120 res*1.4 res])

subplot(234)
plot_sc(U_true_c2,model_true);
% title('True density','interpreter', 'latex');
set(gca,'fontsize',fontsize)
title('True denstiy $\rho_2$','interpreter', 'latex');
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[100  120 res*1.4 res])

% simulation
subplot(232)
plot_sc(U_sim_c1,model_true);
% title('True density','interpreter', 'latex');
set(gca,'fontsize',fontsize)
title('Foward sim $\bar{\rho_1}$','interpreter', 'latex');
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[100  120 res*1.4 res])

subplot(235)
plot_sc(U_sim_c2,model_true);
% title('True density','interpreter', 'latex');
set(gca,'fontsize',fontsize)
title('Foward sim $\bar{\rho_2}$','interpreter', 'latex');
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[100  120 res*1.4 res])


% estimated
subplot(233)
plot_sc(U_est_c1,model_true);
% title('Estimated density','interpreter', 'latex');
set(gca,'fontsize',fontsize)
title('Filter estimate $\hat{\rho_1}$','interpreter', 'latex');
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[100  120 res*1.4 res])

subplot(236)
plot_sc(U_est_c2,model_true);
% title('Estimated density','interpreter', 'latex');
set(gca,'fontsize',fontsize)
title('Filter estimate $\hat{\rho_2}$','interpreter', 'latex');
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[100  120 res*1.4 res])

set(gcf,'position',[100  100 res*2 res*1.1])
end

