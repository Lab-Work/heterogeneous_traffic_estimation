% plots for presentation and paper
%% plot velocity function
close all;figure;
r = [0:0.1:1.8];
v_1 = model_true.vm1-model_true.rm2*r;
v_2 = model_true.vm1-model_true.rm1*r;
plot(r,v_1,'linewidth',2); hold on
plot(r,v_2,'linewidth',2);
ylim([0 2]);
xlim([0 2]);
legend({'V_1(r)','V_2(r)'});
set(gca,'XTick',[0 model_true.rm2 model_true.rm1],'FontName','Times','fontsize',26);
set(gca,'XTickLabel',{'0','r_2^m','r_1^m'},'fontsize',26) 
xlabel('Total density','fontsize',26);
set(gca,'YTick',0:1.8:2,'FontName','Times');
set(gca,'YTickLabel',{'0','v^m'},'fontsize',26) 
ylabel('Velocity','fontsize',26)

%% plot effective particle size histogram
figure;
histogram(N_eff,20,'FaceColor','k','FaceAlpha',0.4)
set(gca,'fontsize',30,'FontName','Times')
xlim([0 800])

%% plot fundamental diagram
figure;
r = [0:0.1:2];
v_1 = model_true.vm1-(model_true.vm1/model_true.rm1)*r;
v_2 = model_true.vm1-(model_true.vm1/model_true.rm2)*r;
q_1 = v_1.*r;
q_2 = v_2.*r;
v_1_app = model_est.vm1-(model_est.vm1/model_est.rm1)*r;
v_2_app = model_est.vm1-(model_est.vm1/model_est.rm2)*r;
q_1_app = v_1_app.*r;
q_2_app = v_2_app.*r;
plot(r,q_1,'linewidth',2,'color','b'); hold on
plot(r,q_2,'linewidth',2,'color','r');
plot(r,q_1_app,'linewidth',2,'color','b','linestyle','-.');
plot(r,q_2_app,'linewidth',2,'color','r','linestyle','-.');
ylim([0 1.2]);
xlim([0 2]);
legend({'True, small','True, large','Approx, small','Approx, large'});
% set(h,'Location','EastOutside',...
%         'position',[0.11 0.80 0.5 0.14],'Orientation','horizontal');
% set(gca,'XTick',[0 model_true.rm2 model_true.rm1],'FontName','Times','fontsize',26);
% set(gca,'XTickLabel',{'0','r_2^m','r_1^m'},'fontsize',26) 
set(gca,'fontsize',26,'FontName','Times') 
xlabel('Total density (r)');
% set(gca,'YTick',0:1.8:2,'FontName','Times');
% set(gca,'YTickLabel',{'0','v^m'},'fontsize',26) 
ylabel('Flow')

%% particle influence
num = xlsread('beeq.xlsx','test3');
close all
figure;
plot(num(2:6,1),num(2:6,2),'k','linewidth',2,'linestyle','-'); % PF
hold on
plot(num(2:6,1),num(10:14,2),'r','linewidth',2,'linestyle','-.'); % SCNM
plot(num(2:6,1),num(18:22,2),'k','linewidth',2,'linestyle','-.'); % PAPF
plot(num(2:6,1),num(26:30,2),'r','linewidth',2,'linestyle','-'); % SCNM+PAPF
xlim([500,2000]); 
xlabel('$N_p$','Interpreter','Latex');
ylabel('Running time (sec/run)')
h = legend({'PF','PF+SCNM','PAPF','PAPF+SCNM'},'Orientation','vertical','fontsize',20);
set(h,'Location','NorthWest');
title('$N_p$ vs. running time','Interpreter','Latex');
set(gca,'fontsize',26,'FontName','Times')
res = 450;
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10  100 res res*0.8])

figure;
plot(num(2:6,1),num(2:6,5),'k','linewidth',2,'linestyle','-'); % PF
hold on
plot(num(2:6,1),num(10:14,5),'r','linewidth',2,'linestyle','-.'); % SCNM
plot(num(2:6,1),num(18:22,5),'k','linewidth',2,'linestyle','-.'); % PAPF
plot(num(2:6,1),num(26:30,5),'r','linewidth',2,'linestyle','-'); % SCNM+PAPF
xlim([500,2000]); 
xlabel('$N_p$','Interpreter','Latex');
ylabel('$N_{\text{eff}}$','Interpreter','Latex');
% h = legend({'PF','SCNM','PAPF','SCNM+PAPF'},'Orientation','vertical','fontsize',20);
set(h,'Location','NorthWest');
title('$N_p$ vs. average $N_{eff}$','Interpreter','Latex');
set(gca,'fontsize',26,'FontName','Times')
res = 450;
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10  100 res res*0.8])

figure;
plot(num(2:6,1),num(2:6,6),'k','linewidth',2,'linestyle','-'); % PF
hold on
plot(num(2:6,1),num(10:14,6),'r','linewidth',2,'linestyle','-.'); % SCNM
plot(num(2:6,1),num(18:22,6),'k','linewidth',2,'linestyle','-.'); % PAPF
plot(num(2:6,1),num(26:30,6),'r','linewidth',2,'linestyle','-'); % SCNM+PAPF
xlim([500,2000]); 
ylim([10 100]);
xlabel('$N_p$','Interpreter','Latex');
ylabel('MAE reduction (\%)','Interpreter','Latex');
% h = legend({'PF','SCNM','PAPF','SCNM+PAPF'},'Orientation','vertical','fontsize',20);set(h,'Location','NorthWest');
title('$N_p$ vs. $\rho_1$ accuracy','Interpreter','Latex');
set(gca,'fontsize',26,'FontName','Times')
res = 450;
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10  100 res res*0.8])

figure;
plot(num(2:6,1),num(2:6,7),'k','linewidth',2,'linestyle','-'); % PF
hold on
plot(num(2:6,1),num(10:14,7),'r','linewidth',2,'linestyle','-.'); % SCNM
plot(num(2:6,1),num(18:22,7),'k','linewidth',2,'linestyle','-.'); % PAPF
plot(num(2:6,1),num(26:30,7),'r','linewidth',2,'linestyle','-'); % SCNM+PAPFxlim([min(num(:,1)) max(num(:,1))]); 
xlim([500,2000]); 
ylim([10 100]);
xlabel('$N_p$','Interpreter','Latex');
ylabel('MAE reduction (\%)','Interpreter','Latex');
% h = legend({'PF','SCNM','PAPF','SCNM+PAPF'},'Orientation','vertical','fontsize',20);set(h,'Location','NorthWest');
title('$N_p$ vs. $\rho_2$ accuracy','Interpreter','Latex');
set(gca,'fontsize',26,'FontName','Times')
res = 450;
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10  100 res res*0.8])

%% plot boundaries
figure;
plot(model_est.d1l);
hold on
plot(model_true.d1l);
legend({'noisy upstream','true upstream'},'Location','southeast');
xlabel('Timesteps');
ylabel('Density');
set(gca,'Fontsize',26);

%% plot initial
figure;
plot(U_est_c2(1,:));
hold on
plot(U_true_c2(1,:));
legend({'approx upstream','true upstream'});
xlabel('# Cell');
ylabel('Density');
set(gca,'Fontsize',26);
%% plot white noise vs. correlated noise
close all
N = 80;
tau = 1:N;
tau_m = triu(toeplitz(tau));
tau_m = tau_m-tau_m';
R_m = autocorrelate(tau_m.^2,20^2);
[V,D] = eig(R_m);
x_miu = zeros(size(tau));
sum_init = 0;
for i = 1:size(V,1)
    sum_init = sum_init + sqrt(D(i,i)) * randn * V(:,i)';
end
corr_noise = sum_init;
white_noise = normrnd(0,1,[1 N]);
figure;
subplot(211)
plot(white_noise);
title('Standard Gaussian white noise','Fontsize',20,'interpreter','latex');
subplot(212)
plot(corr_noise);
title('Standard Gaussian random process','Fontsize',20,'interpreter','latex');

%% plot autocorrelation function
tau = -600:1:600;
R = autocorrelate(tau,100);
plot(tau,R,'linewidth',2);
xlabel('$\tau$','Fontsize',20,'interpreter','latex');
ylabel('$R(\tau)$','Fontsize',20,'interpreter','latex');
title('Autocorrelation function','Fontsize',20,'interpreter','latex');

%% plot group stacked bars
% clc;
% close all
% figure
% subplot(121)
% stackData(1,:,:) = [0.1086,0.038;0.0762,0.0584];
% stackData(2,:,:) = [0.0448,0.034;0.0607,0.0511];
% groupLabels = {'LWR model','Creeping model'};
% barLabels = {'',''};
% plotBarStackGroups(stackData, groupLabels,barLabels);
% title('Creeping test','fontsize',20);
% 
% subplot(122)
% tackData(1,:,:) = [0.0819,0.0345;0.0815,0.0525];
% stackData(2,:,:) = [0.0572,0.0311;0.0562,0.0546];
% groupLabels = {'LWR model','Creeping model'};
% barLabels = {'',''};
% plotBarStackGroups(stackData, groupLabels,barLabels);
% title('Overtaking test','fontsize',20);
clc;
close all
figure
subplot(121)
stackData(1,:,:) = [0.1086,0.038;0.0762,0.0584];
stackData(2,:,:) = [0.0448,0.034;0.0607,0.0511];
groupLabels = {'LWR model','Creeping model'};
barLabels = {'',''};
plotBarStackGroups(stackData, groupLabels,barLabels);
title('Creeping test','fontsize',20);

subplot(122)
tackData(1,:,:) = [0.0819,0.0345;0.0815,0.0525];
stackData(2,:,:) = [0.0572,0.0311;0.0562,0.0546];
groupLabels = {'LWR model','Creeping model'};
barLabels = {'',''};
plotBarStackGroups(stackData, groupLabels,barLabels);
title('Overtaking test','fontsize',20);
%% plot group stacked bars IEEE conf
% clc;
% close all
% figure
% stackData(1,:,:) = [0.0582,0.0253;0.0733,0.0696];
% stackData(2,:,:) = [0.0627,0.0229;0.0679,0.0281];
% groupLabels = {'Case 1','Case 2'};
% barLabels = {'',''};
% plotBarStackGroups(stackData, groupLabels,barLabels);
% % title('Creeping test','fontsize',20);

clc;
close all
figure
stackData(1,:,:) = [0.081,0.012;0.106,0.053];
stackData(2,:,:) = [0.048,0.027;0.051,0.047];
groupLabels = {'Creeping','Overtaking'};
barLabels = {'',''};
plotBarStackGroups(stackData, groupLabels,barLabels);
% title('Creeping test','fontsize',20);

%% sensitivity analysis
[num,txt,raw] = xlsread('../../beeq.xlsx', 4,'A80:H111');
% num = xlsread('error.xlsx', 6);
test1 = 1:5;
test2 = 9:13;
test3 = 17:21;
test4 = 25:29;

% baseline values
val1 = [32.88 31.57];
val2 = [18.45 61.57];
val3 = [50.49 81.39];
val4 = [15.77 41.84];

xa = test3;
true = val3;

close all
figure;
subplot(211)
plot(num(xa,1),num(xa,2),'-k','linewidth',2); % vm1
hold on
plot(num(xa,1),num(xa,3),':k','linewidth',2);
plot(num(xa,1),num(xa,4),'-.k','linewidth',2);
plot(num(xa,1),true(1)*ones(1,5),'-r','linewidth',2);
ylim([-5 70]);
% set(gca,'XTick',0:1.8:2,'FontName','Times');
set(gca,'XTickLabel','','fontsize',26,'FontName','Times')
% xlabel('Perturbation','fontsize',26);
% set(gca,'YTick',0:1.8:2,'FontName','Times');
% set(gca,'YTickLabel',{'0','v^m'},'fontsize',20) 
ylabel('$\rho_1$ ','interpreter','Latex')
% h = legend({'v^m','r_1^m','r_2^m','no perturbation'},'Orientation','horizontal','fontsize',20);
% set(h,'Location','NorthWest',...
%     'position',[0.38 0.79 0.3 0.08]);
title('Percentage of MAE reduction (\%)','interpreter','latex');

subplot(212)
plot(num(xa,1),num(xa,6),'-k','linewidth',2); hold on
plot(num(xa,1),num(xa,7),':k','linewidth',2);
plot(num(xa,1),num(xa,8),'-.k','linewidth',2);
plot(num(xa,1),true(2)*ones(1,5),'-r','linewidth',2);

ylim([20 90]);
% set(gca,'XTick',0:1.8:2,'FontName','Times');
set(gca,'XTickLabel',{'-10%','-5%','0','5%','10%'},'fontsize',26,'FontName','Times') 
xlabel('Perturbation','fontsize',26);
% set(gca,'YTick',0:1.8:2,'FontName','Times');
% set(gca,'YTickLabel',{'0','v^m'},'fontsize',20) 
ylabel('$\rho_2$','interpreter','Latex')
% h = legend({'v^m','r_1^m','r_2^m'},'Orientation','horizontal','fontsize',20);
% set(h,'Location','NorthWest',...
%     'position',[0.38 0.32 0.3 0.08])
% title('Percentage of MAE reduction (\%)','interpreter','latex');
set(gcf,'position',[500 300 450 450])

%% influence of particles
% save runtime
for n = 1:100
    A = rand(n,n);
    b = rand(n,1);
    tic
    x = A\b;
    t(n) = toc;
end
plot(t)

%% plot MAE greyscale
figure;
subplot(221)
imagesc(abs(U_true_c1-U_sim_c1));
colormap(gray)
colorbar
caxis([0 0.7]);
set(gca,'YDir','normal');

subplot(222)
imagesc(abs(U_true_c1-U_est_c1));
colormap(gray)
colorbar
caxis([0 0.7]);
set(gca,'YDir','normal');

subplot(223)
imagesc(abs(U_true_c2-U_sim_c2));
colormap(gray)
colorbar
caxis([0 0.7]);
set(gca,'YDir','normal');

subplot(224)
imagesc(abs(U_true_c2-U_est_c2));
colormap(gray)
colorbar
caxis([0 0.7]);
set(gca,'YDir','normal');

plot_sc_gray()