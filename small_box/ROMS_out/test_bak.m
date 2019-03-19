
clear
clc
proj_root = '~/Documents/GitHub/Rutgers_ROMS/projects_ROMS/small_box';

%%

time   = ncread('roms_his.nc','ocean_time');
dstart = ncread('roms_his.nc','dstart'); % days since t_ref

zeta    = ncread('roms_his.nc','zeta');
ubar    = ncread('roms_his.nc','ubar');
vbar    = ncread('roms_his.nc','vbar');
u       = ncread('roms_his.nc','u');
v       = ncread('roms_his.nc','v');
temp    = ncread('roms_his.nc','temp');
salt    = ncread('roms_his.nc','salt');

% background vertical mixing coefficient for tracers (temp & salt)
Akt_bak = ncread('roms_his.nc','Akt_bak');
AKs     = ncread('roms_his.nc','AKs'); % vertical mxing coefficient for salt
AKt     = ncread('roms_his.nc','AKt'); % vertical mxing coefficient for temp
AKv     = ncread('roms_his.nc','AKv'); % vertical mxing coefficient for momentum
tke     = ncread('roms_his.nc','tke');
gls     = ncread('roms_his.nc','gls'); % turbulent generic length scale

t_ref = datenum('0001-01-01 00:00:00','yyyy-mm-dd HH:MM:SS');
time  = time/3600/24 + t_ref;

%%

Vtransform  = ncread('roms_his.nc','Vtransform');
Vstretching = ncread('roms_his.nc','Vstretching');
theta_s     = ncread('roms_his.nc','theta_s');
theta_b     = ncread('roms_his.nc','theta_b');
hc          = ncread('roms_his.nc','hc');
h           = ncread('roms_his.nc','h');
%grid        = ncread('roms_his.nc','grid');

N     = 180;
igrid = 1; % for RHO points

z_rho = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
                  igrid, h, zeta(:,:,1));

%% Extraction

temp_ini = temp(:,:,:,1);
temp_end = temp(:,:,:,end);

salt_ini = salt(:,:,:,1);
salt_end = salt(:,:,:,end);

ubar_ini = ubar(:,:,1);
ubar_end = ubar(:,:,end);

vbar_ini = vbar(:,:,1);
vbar_end = vbar(:,:,end);

u_Hov = squeeze(squeeze(u(1,1,:,:)));
v_Hov = squeeze(squeeze(v(1,1,:,:)));
temp_Hov = squeeze(squeeze(temp(1,1,:,:)));
salt_Hov = squeeze(squeeze(salt(1,1,:,:)));

AKt_Hov = squeeze(squeeze(AKt(1,1,:,:)));
AKv_Hov = squeeze(squeeze(AKv(1,1,:,:)));
tke_Hov = squeeze(squeeze(tke(1,1,:,:)));
gls_Hov = squeeze(squeeze(gls(1,1,:,:)));

%%

figure('position', [0, 0, 500, 600]);
plot(squeeze(temp_ini(1,1,:)),squeeze(z_rho(1,1,:)))
hold on
plot(squeeze(temp_end(1,1,:)),squeeze(z_rho(1,1,:)))
hold off
legend({'temp. - initial','temp. - after 15 months'},'Location','east',...
    'FontSize',11,'Interpreter','latex')
ylim([-300 0])

% set(gca,'LooseInset', get(gca,'TightInset')); % no blank edge
% saveas(gcf, [proj_root,'/Figs/tprof'], 'png');

%%

figure('position', [0, 0, 500, 600]);
plot(squeeze(salt_ini(1,1,:)),squeeze(z_rho(1,1,:)))
hold on
plot(squeeze(salt_end(1,1,:)),squeeze(z_rho(1,1,:)))
hold off
legend({'sal. - initial','sal. - after 15 months'},'Location','best',...
    'FontSize',11,'Interpreter','latex')
ylim([-300 0])

% set(gca,'LooseInset', get(gca,'TightInset')); % no blank edge
% saveas(gcf, [proj_root,'/Figs/sprof'], 'png');
