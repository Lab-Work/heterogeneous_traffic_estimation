function plot_six_gray(U_true_c1,U_sim_c1,U_est_c1,U_true_c2,U_sim_c2,U_est_c2,model_true)
res = 600;
fontsize = 20;

figure;
subplot(231)
plot_sc(U_true_c1,model_true);
title('True denstiy $\rho_1$','interpreter', 'latex');
colorbar;
set(gca,'fontsize',fontsize)
set(gca,'xlabel',[])

subplot(234)
plot_sc(U_true_c2,model_true);
set(gca,'fontsize',fontsize)
title('True denstiy $\rho_2$','interpreter', 'latex');


% simulation
subplot(232)
plot_sc(abs(U_true_c1-U_sim_c1),model_true);
colormap(gray)
set(gca,'YDir','normal');
set(gca,'fontsize',fontsize)
set(gca,'xlabel',[])
set(gca,'ylabel',[])
title('$\bar{e}_1$','interpreter', 'latex');

subplot(235)
plot_sc(abs(U_true_c2-U_sim_c2),model_true);
set(gca,'fontsize',fontsize)
set(gca,'ylabel',[])
title('$\bar{e}_2$','interpreter', 'latex');


% estimated
subplot(233)
plot_sc(abs(U_true_c1-U_est_c1),model_true);
set(gca,'fontsize',fontsize)
set(gca,'xlabel',[])
set(gca,'ylabel',[])
title('$\hat{e}_1$','interpreter', 'latex');

subplot(236)
plot_sc(abs(U_true_c2-U_est_c2),model_true);
set(gca,'fontsize',fontsize)
title('$\hat{e}_2$','interpreter', 'latex');
set(gca,'ylabel',[])
set(gcf,'position',[100  100 res*1.6 res*0.8])
end

