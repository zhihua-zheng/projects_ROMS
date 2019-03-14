
clear

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
tke     = ncread('roms_his.nc','tke');
gls     = ncread('roms_his.nc','gls');
AKs     = ncread('roms_his.nc','AKs');
AKt     = ncread('roms_his.nc','AKt');
Akt_bak = ncread('roms_his.nc','Akt_bak');

t_ref = datenum('0001-01-01 00:00:00','yyyy-mm-dd HH:MM:SS');
time  = time/3600/24 + t_ref;

%%

Vtransform  = ncread('roms_his.nc','Vtransform');
Vstretching = ncread('roms_his.nc','Vstretching');
theta_s     = ncread('roms_his.nc','theta_s');
theta_b     = ncread('roms_his.nc','theta_b');
hc          = ncread('roms_his.nc','hc');
h           = ncread('roms_his.nc','h');
grid        = ncread('roms_his.nc','grid');

N     = 180;
igrid = 1; % for RHO points

z_rho = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
                  igrid, h, zeta);


%%

temp_ini = temp(:,:,:,1);
temp_end = temp(:,:,:,end);

plot(squeeze(temp_ini(1,1,:)),z_rho);
