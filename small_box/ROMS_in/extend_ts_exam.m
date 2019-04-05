%% extend_ts_exam

figure('position', [0, 0, 400, 700])
scatter(t_ic,z_ic,25,'filled'); axis ij
hold on 
plot(t_new,z_new,'Linewidth',2)
hold off 
box on
export_fig('./figs/temp','-png','-transparent','-painters')

figure('position', [0, 0, 400, 700])
scatter(s_ic,z_ic,15,'filled'); axis ij
hold on 
plot(s_new,z_new,'Linewidth',2)
hold off 
box on
export_fig('./figs/sal','-png','-transparent','-painters')

%----- Gradient
t_new_z = -center_diff(t_new,z_new,1); % dT/dz
s_new_z = -center_diff(s_new,z_new,1); % dS/dz
z_mid_new = (z_new(1:end-1) + z_new(2:end))/2;

t_ic_z = -center_diff(t_ic,z_ic,1); % dT/dz in original profile
s_ic_z = -center_diff(s_ic,z_ic,1); % dS/dz in original profile
z_mid_ic = (z_ic(1:end-1) + z_ic(2:end))/2;

figure('position', [0, 0, 400, 700])
line(s_new_z,z_mid_new,'Linewidth',1.8,'Color',rgb('vermillion')) 
line(t_new_z,z_mid_new,'Linewidth',1.8,'Color',rgb('soft blue'))
line([0 0],[0 4500],'Color',[.3 .3 .3],'LineStyle','--','Linewidth',1)
hold on
scatter(s_ic_z,z_mid_ic,7,'filled','ks')
hold on
scatter(t_ic_z,z_mid_ic,7,'filled','ks')
axis ij
box on
legend({'$dS/dz$','$dT/dz$'},'Interpreter','latex','fontsize',20,...
    'Location','best','color','none')
export_fig('./figs/ts_gradient','-png','-transparent','-painters')
%----- Gradient
