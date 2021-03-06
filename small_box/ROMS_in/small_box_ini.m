%% small_box_init

% Script to create initial conditions for the ROMS small box test case
% Zhihua Zheng, UW-APL, Mar. 5 2019

clear

%% 

show_fig = 0; % turn on/off figure section

%  Set input/output NetCDF files.
my_root = '~/Documents/GitHub/Rutgers_ROMS/projects_ROMS/small_box';

GRDname      = fullfile(my_root,'Data/grid','uniform_grid_50m.nc');
CLname       = fullfile(my_root,'Data/CL',  'ocsp_cl_1995_2017.mat');
LES_INIname  = fullfile(my_root,'Data/SS1D/nps10f192','nps10f192MYinit.dat');
LES_STKname  = fullfile(my_root,'Data/SS1D/nps10f192','nps10f192_my_stksdf.dat');
LES_SDIRname = fullfile(my_root,'Data/SS1D/nps10f192','nps10f192_my_sdirad.dat');
% INIname      = fullfile(my_root,'ROMS_in','small_box_ini.nc');
INIname     = 'small_box_ini.nc';

%% Box configuration

%--------------------------------------------------------------------------
%  Set application parameters in structure array, S.
%--------------------------------------------------------------------------

S.ncname        = INIname; % initial conditions output file name.
S.spherical     = 0; % use Cartesian grid

% for RHO-point
S.Lm            = 3; % number of interior RHO-points in X-direction
S.Mm            = 3; % number of interior RHO-points in Y-direction
S.N             = 180; % number of vertical levels at RHO-points

S.NT            = 2; % number of tracers

% vertical S-coordinate and vertical stretching function
S.Vtransform    = 2;
S.Vstretching   = 4; % (Shchepetkin 2010)
S.theta_s       = 10; 
S.theta_b       = 0; 
S.Tcline        = 240; 
S.hc            = S.Tcline;

%--------------------------------------------------------------------------
%  Create initial condition Netcdf file.
%--------------------------------------------------------------------------

[~]=c_initial(S);

%  Set attributes for "ocean_time".

avalue='seconds since 0001-01-01 00:00:00';
[~]=nc_attadd(INIname,'units',avalue,'ocean_time');
  
avalue='360.0 days in every year';
[~]=nc_attadd(INIname,'calendar',avalue,'ocean_time');

%--------------------------------------------------------------------------
%  Set grid variables.
%--------------------------------------------------------------------------

V=nc_vnames(GRDname);
nvars=length(V.Variables);

%  Horizontal grid variables. Read in for input GRID NetCDF file.
if (S.spherical)
  S.lon_rho = nc_read(GRDname, 'lon_rho');
  S.lat_rho = nc_read(GRDname, 'lat_rho');
  
  S.lon_u   = nc_read(GRDname, 'lon_u');
  S.lat_u   = nc_read(GRDname, 'lat_u');
  
  S.lon_v   = nc_read(GRDname, 'lon_v');
  S.lat_v   = nc_read(GRDname, 'lat_v');
else  
  S.x_rho   = nc_read(GRDname, 'x_rho');
  S.y_rho   = nc_read(GRDname, 'y_rho');
  
  S.x_u     = nc_read(GRDname, 'x_u');
  S.y_u     = nc_read(GRDname, 'y_u');
  
  S.x_v     = nc_read(GRDname, 'x_v');
  S.y_v     = nc_read(GRDname, 'y_v');  
end  

%  Read in Land/Sea mask, if appropriate.
for n=1:nvars
  name=char(V.Variables(n).Name);
  switch (name)
    case 'mask_rho'
      S.mask_rho = nc_read(GRDname, 'mask_rho');
    case 'mask_u'
      S.mask_u   = nc_read(GRDname, 'mask_u');
    case 'mask_v'
      S.mask_v   = nc_read(GRDname, 'mask_v');
  end
end

%  Bathymetry
S.h = nc_read(GRDname, 'h');

%  Set vertical grid variables.
[S.s_rho,S.Cs_r] = stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N, 0, 0);
[S.s_w,S.Cs_w]   = stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N, 1, 0);
                        
%--------------------------------------------------------------------------
%  Set zero initial conditions.
%--------------------------------------------------------------------------

Lr = S.Lm+2;   Lu = Lr-1;   Lv = Lr; % X-direction [\xi]
Mr = S.Mm+2;   Mu = Mr;     Mv = Mr-1; % Y-direction [\eta]

S.zeta = zeros([Lr Mr]);
S.ubar = zeros([Lu Mu]);
S.vbar = zeros([Lv Mv]);
S.u    = zeros([Lu Mu S.N]);
S.v    = zeros([Lv Mv S.N]);
S.temp = zeros([Lr Mr S.N]);
S.salt = zeros([Lr Mr S.N]);

%  If Land/Sea masking arrays are not found, initialize them to unity.
if (~isfield(S, 'mask_rho')),  S.mask_rho = ones([Lr Mr]);  end,
if (~isfield(S, 'mask_u'  )),  S.mask_u   = ones([Lu Mu]);  end,
if (~isfield(S, 'mask_v'  )),  S.mask_v   = ones([Lv Mv]);  end,

%--------------------------------------------------------------------------
%  Write out grid variables.
%--------------------------------------------------------------------------
			 
[~]=nc_write(INIname,   'spherical',   S.spherical);

[~]=nc_write(INIname,   'Vtransform',  S.Vtransform);
[~]=nc_write(INIname,   'Vstretching', S.Vstretching);
[~]=nc_write(INIname,   'theta_s',     S.theta_s);
[~]=nc_write(INIname,   'theta_b',     S.theta_b);
[~]=nc_write(INIname,   'Tcline',      S.Tcline);
[~]=nc_write(INIname,   'hc',          S.hc);

[~]=nc_write(INIname,   's_rho',       S.s_rho);
[~]=nc_write(INIname,   's_w',         S.s_w);
[~]=nc_write(INIname,   'Cs_r',        S.Cs_r);
[~]=nc_write(INIname,   'Cs_w',        S.Cs_w);

[~]=nc_write(INIname,   'h',           S.h);

if (S.spherical)
  [~]=nc_write(INIname, 'lon_rho',     S.lon_rho);
  [~]=nc_write(INIname, 'lat_rho',     S.lat_rho);
  [~]=nc_write(INIname, 'lon_u',       S.lon_u);
  [~]=nc_write(INIname, 'lat_u',       S.lat_u);
  [~]=nc_write(INIname, 'lon_v',       S.lon_v);
  [~]=nc_write(INIname, 'lat_v',       S.lat_v);
else
  [~]=nc_write(INIname, 'x_rho',       S.x_rho);
  [~]=nc_write(INIname, 'y_rho',       S.y_rho);
  [~]=nc_write(INIname, 'x_u',         S.x_u);
  [~]=nc_write(INIname, 'y_u',         S.y_u);
  [~]=nc_write(INIname, 'x_v',         S.x_v);
  [~]=nc_write(INIname, 'y_v',         S.y_v);
end

%--------------------------------------------------------------------------
%  Compute depths at horizontal and vertical RHO-points.
%--------------------------------------------------------------------------

igrid = 1; % for RHO points
z_rho = set_depth(S.Vtransform, S.Vstretching,                          ...
                  S.theta_s, S.theta_b, S.hc, S.N,                      ...
                  igrid, S.h, S.zeta);

%% Read in LES initial profiles and Climatology

load(CLname)
icdata = load(LES_INIname);

% Note the LES initial velocity includes Stokes drift, the coordinate
% system is directly reversed from atmospheric boundary layer LES model,
% and LES wind is imposed towards North

z_ic = icdata(:,1);
t_ic = icdata(:,2);
s_ic = icdata(:,3);
u_ic = icdata(:,4);  % [m/s]
v_ic = icdata(:,5);
% q2_ic  = icdata(:,6);
% q2l_ic = icdata(:,7);

n      = length(z_ic);
dz_LES = 1.42; % see HD2008

%% Stokes drift profile in LES intial condition

st_data = load(LES_STKname);
st_dir  = load(LES_SDIRname);

f_ctr  = st_data(1,2:end); % frequencies
stksdf = st_data(2,2:end); % Stokes spectrum, dU_st
sdirad = st_dir(2,2:end); % Stokes drift direction

% the x, y Stokes components in LES
du_st = stksdf .* sin(sdirad);
dv_st = stksdf .* cos(sdirad);

% get the Stokes drift profile
g         = 9.81;
k_ctr     = (2*pi*f_ctr).^2/g;

% multiplier produced by spatial filter of exp(2kz)
filter_m  = sinh(k_ctr*dz_LES) ./ (k_ctr*dz_LES);

z_decay_u = exp(-2*z_ic*k_ctr);
z_decay_v = exp(-2*z_ic*k_ctr);
u_st      = z_decay_u*(du_st .* filter_m)';
v_st      = z_decay_v*(dv_st .* filter_m)';
% Wave's direction are in usual coordinate

% velocity conversion, the wind is imposed towards East in ROMS
%   u_roms = v_LES = - v_ic - v_st
% - v_roms = u_LES =   u_ic - u_st

% subtract the Stokes drift from total velocity & adjust y-axis in LES
u_LES  =  u_ic - u_st;
v_LES  = -v_ic - v_st;

u_roms =  v_LES;
v_roms = -u_LES;

%% Extend LES initial profiles into depth

temp_cl_z = center_diff(temp_cl,z,2); % gradient towards bottom
sal_cl_z  = center_diff(sal_cl,z,2); % gradient towards bottom
z_mid = (z(1:end-1) + z(2:end))/2;

z_cutoff = z_ic(end);
[~,i_start] = find(z_mid-z_cutoff > 0,1,'first');

t_z_cutoff = (t_ic(end)-t_ic(end-1))/(z_cutoff-z_ic(end-1));
s_z_cutoff = (s_ic(end)-s_ic(end-1))/(z_cutoff-z_ic(end-1));

t_z = [ones(1,n-1)*NaN t_z_cutoff temp_cl_z(i_start:end)]';
s_z = [ones(1,n-1)*NaN s_z_cutoff sal_cl_z(i_start:end)]';

z_new = [z_ic; z_mid(i_start:end)'];
t_new = zeros(size(z_new));
s_new = zeros(size(z_new));
u_new = zeros(size(z_new));
v_new = zeros(size(z_new));

u_new(1:n) = u_roms;
v_new(1:n) = v_roms;
t_new(1:n) = t_ic;
s_new(1:n) = s_ic;

t_deep = cumtrapz(z_new(n:end), t_z(n:end));
s_deep = cumtrapz(z_new(n:end), s_z(n:end));

t_new(n:end) = t_ic(end) + t_deep;
s_new(n:end) = s_ic(end) + s_deep;

%% Display

if show_fig
    extend_ts_exam;
end

%% Interpolate profiles to ROMS-z level

for i = 1:Lr
    for j = 1:Mr
        
        S.temp(i,j,:) = interp1(-z_new,t_new,z_rho(i,j,:));
        S.salt(i,j,:) = interp1(-z_new,s_new,z_rho(i,j,:));   
    end
end

for i = 1:Lu
    for j = 1:Mu
        
        S.ubar(i,j)   = trapz(-z_new,u_new) / S.h(i,j);
        S.u(i,j,:)    = interp1(-z_new,u_new,z_rho(i,j,:));
    end
end

for i = 1:Lv
    for j = 1:Mv
        
        S.vbar(i,j)   = trapz(-z_new,v_new) / S.h(i,j);
        S.v(i,j,:)    = interp1(-z_new,v_new,z_rho(i,j,:));
    end
end

%%  Write out initial conditions

IniRec = 1;                               % NetCDF time record - Tindex

S.ocean_time = 30.0*86400;                % initial conditions time (s)

[~]=nc_write(INIname, 'ocean_time', S.ocean_time, IniRec);

[~]=nc_write(INIname, 'zeta', S.zeta, IniRec);
[~]=nc_write(INIname, 'ubar', S.ubar, IniRec);
[~]=nc_write(INIname, 'vbar', S.vbar, IniRec);
[~]=nc_write(INIname, 'u',    S.u,    IniRec);
[~]=nc_write(INIname, 'v',    S.v,    IniRec);
[~]=nc_write(INIname, 'temp', S.temp, IniRec);
[~]=nc_write(INIname, 'salt', S.salt, IniRec);
