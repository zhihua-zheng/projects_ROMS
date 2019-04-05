

%% General

clear
clc

if ismac
    proj_root = '/Volumes/projects_ROMS/small_box';
elseif isunix
    proj_root = '~/Documents/GitHub/Rutgers_ROMS/projects_ROMS/small_box';
% elseif ispc
    % Code to run on Windows platform
else
    disp('Platform not supported')
end

%HISname = [proj_root,'/ROMS_out/Inertial/roms_his.nc'];
HISname = [proj_root,'/ROMS_out/TransWind/roms_his.nc'];

%%

time   = ncread(HISname,'ocean_time');
dt     = ncread(HISname,'dt');
nHIS   = ncread(HISname,'nHIS');
dstart = ncread(HISname,'dstart'); % days since t_ref

sustr = ncread(HISname,'sustr'); 
svstr = ncread(HISname,'svstr');

zeta    = ncread(HISname,'zeta');
ubar    = ncread(HISname,'ubar');
vbar    = ncread(HISname,'vbar');
u       = ncread(HISname,'u');
v       = ncread(HISname,'v');
temp    = ncread(HISname,'temp');
salt    = ncread(HISname,'salt');

% background vertical mixing coefficient for tracers (temp & salt)
Akt_bak = ncread(HISname,'Akt_bak');
AKs     = ncread(HISname,'AKs'); % vertical mxing coefficient for salt
AKt     = ncread(HISname,'AKt'); % vertical mxing coefficient for temp
AKv     = ncread(HISname,'AKv'); % vertical mxing coefficient for momentum
tke     = ncread(HISname,'tke');
gls     = ncread(HISname,'gls'); % turbulent generic length scale

t_ref = datenum('0001-01-01 00:00:00','yyyy-mm-dd HH:MM:SS');
time  = time/3600/24 + t_ref;


%% Wind check

figure('position', [0, 0, 600, 120])
plot(time,squeeze(sustr(3,3,:))); datetick('x','dd'); xlim(time([1,end]))
hold on
plot(time,squeeze(svstr(3,3,:))); datetick('x','dd'); xlim(time([1,end]))

%% Grid

Vtransform  = ncread(HISname,'Vtransform');
Vstretching = ncread(HISname,'Vstretching');
theta_s     = ncread(HISname,'theta_s');
theta_b     = ncread(HISname,'theta_b');
hc          = ncread(HISname,'hc');
h           = ncread(HISname,'h');
f_Coriolis  = ncread(HISname,'f');

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

ubar_a = squeeze(ubar(3,3,:));
vbar_a = squeeze(vbar(3,3,:));

u_surf = squeeze(u(3,3,end,:));
v_surf = squeeze(v(3,3,end,:));

u_Hov = squeeze(squeeze(u(3,3,:,:)));
v_Hov = squeeze(squeeze(v(3,3,:,:)));
temp_Hov = squeeze(squeeze(temp(3,3,:,:)));
salt_Hov = squeeze(squeeze(salt(3,3,:,:)));

AKt_Hov = squeeze(squeeze(AKt(3,3,:,:)));
AKv_Hov = squeeze(squeeze(AKv(3,3,:,:)));
tke_Hov = squeeze(squeeze(tke(3,3,:,:)));
gls_Hov = squeeze(squeeze(gls(3,3,:,:)));

%% T-profile

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

%% S-profile

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

%% Velocity

cur_a = complex(ubar_a,vbar_a);

t_Coriolis = 2*pi/f_Coriolis(1)/3600; % [hour]

%------- Welch's power spectral density estimate --------------------------

cur_a = cur_a - mean(cur_a);
 
% the next power of 2 greater than signal length - FFT transfrom length
n = 2^nextpow2(length(cur_a));

% the number of samples per unit time [day]
fs = 24*3600/(nHIS*dt); 

[p_cur, f] = pwelch(cur_a,[],[],n,fs); % returned f in cycle/day

rotary_spec(f,p_cur,24/t_Coriolis,0)
%--------------------------------------------------------------------------

