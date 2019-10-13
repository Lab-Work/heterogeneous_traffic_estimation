function plot_sc(rho,model)
% plot the heat map of density over space and time
h = imagesc(rho(:, 2:end-1));
set(gca,'YDir','normal');
%     xticklabels([10:10:model.N-1]);
a = get(gca,'XTick');
set(gca,'XTick',a,'FontName','Times')

colorbar;
caxis([0 1.0])
xlabel(sprintf('Space $\\Delta x = %.2f$', model.dx),'interpreter', 'latex','Fontsize',20);
ylabel(sprintf('Time $\\Delta t = %.2f$', model.dt),'Interpreter','latex')
% pos = get(gca, 'Position');
% set(gca, 'Position', [pos(1)+0.04 pos(2)+0.08 pos(3) pos(4)-0.08]);
end

