function plot_eps(N_eff)
% plot the effective particle size in histogram

figure;
[counts,centers] = hist(N_eff);
bar(centers,counts,'facecolor','k','edgecolor','w','facealpha',0.5);
xlabel('Effective particle size','interpreter', 'latex'); ylabel('Number of timesteps','interpreter', 'latex');
set(gca,'FontSize',20)
a = get(gca,'XTick');
set(gca,'XTick',a,'FontName','Times');

end

