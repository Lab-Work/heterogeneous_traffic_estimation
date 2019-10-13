function fig = plot_est(n,U_true,U_est,model,pf,U_meas_true,x_next)
% plot the true state and the estimated state
fig = subplot(2,1,1);
%     plot(model.x(2:1:end),U_true{1}(1,2:1:end),'--','color',[0.8,0,0],'linewidth',2), hold on
plot(model.x(2:1:end),U_true{n}(1,2:end),'-','color',[.6,.6,.6],'linewidth',3), hold on
plot(model.x(2:end),U_est{n}(1,2:end),'-','color',[.5,.7,.3],'linewidth',3)
plot(model.x(pf.meas_pt),U_meas_true{n}(1,:),'r*','markers',20);
for i = 1:model.N
    if mod(i,model.N/20)==0
        scatter(model.x(i).*ones(pf.Np,1),x_next(1,i,:),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.1)
        hold on
    end
end
hold off

%     title(sprintf('Small vehicle density at t=%d',n),'interpreter', 'latex')
title(sprintf('timestep=%d',n),'interpreter', 'latex');
axis([0 model.len 0 1.5])
%     h = legend('Initial state','True density','Estimated state','Measurement data');
% h = legend('True density','Estimated state','Sensor location','Particle estimates');
%     set(h,'Location','NorthWest','interpreter', 'latex',...
%         'position',[0.11 0.80 0.5 0.14]);
%     set(h,'Location','NorthWest','interpreter', 'latex','fontsize',22);
ylim([0,2]);
set(gca,'linewidth',1.2,'FontName','Times')
set(gca,'fontsize',30)
set(gca,'xtick',[])
%     set(gca,'position',[.09 .54 .87 .41])
set(gca,'position',[.14 .53 .84 .38])
ylabel('$\rho_1$','interpreter', 'latex')

subplot(2,1,2),
%     plot(model.x(2:1:end),U_true{1}(2,2:1:end),'--','color',[0.8,0,0],'linewidth',2), hold on % initial
plot(model.x(2:end),U_true{n}(2,2:end),'-','color',[0.6,0.6,0.6],'linewidth',3), hold on % true state
plot(model.x(2:end),U_est{n}(2,2:end),'-','color',[.5,.7,.3],'linewidth',3) % estimated state

plot(model.x(pf.meas_pt),U_meas_true{n}(2,:),'r*','markers',20);
for i = 1:model.N
    if mod(i,model.N/20)==0
        scatter(model.x(i).*ones(pf.Np,1),x_next(2,i,:),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',.1)
        hold on
    end
end
hold off

%     title(sprintf('Large vehicle density at t=%d',n),'interpreter', 'latex')
%     h = legend('Initial state','True density','Estimated state','Measurement data');
% h = legend('True density','Estimated state','Sensor location','Particle estimates');
%     set(h,'Location','NorthWest','interpreter', 'latex',...
%         'position',[0.11 0.33 0.5 0.14])
set(gca,'linewidth',1.2)
set(gca,'fontsize',30)
axis([0 model.len 0 1.4])

ylabel('$\rho_2$','interpreter', 'latex')
xlabel('space','interpreter', 'latex');
xh = get(gca,'xlabel'); % handle to the label object
p = get(xh,'position'); % get the current position property
p(2) = p(2)*0.7;       % negative values put the label below the axis
set(xh,'position',p)   % set the new position
%     set(gca,'position',[.09 .07 .87 .41])
set(gca,'position',[.14 .15 .84 .30],'FontName','Times')
%     res = 800;
%     set(gcf,'paperpositionmode','auto')
%     set(gcf,'position',[10  50 res res*1.1])
res = 600;
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10  100 res res*0.70])

end