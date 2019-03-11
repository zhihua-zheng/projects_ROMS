%This is the master code to create an intial linear stratification profile
%followed by creation of a the roms initial file

addpath(genpath('/Users/n2kumar/MAT_MyToolbox/MATLAB_ROMS/'));
sname   = 'initial_temperature_salinity_velocity_zeta_hbeach.mat';

%% STEP 1: Get the NETCDF grid
fname   = '/Users/n2kumar/Documents/MM_Project/Bathymetry/hbeach_grid.nc';
x_rho   = ncread(fname,'x_rho');
y_rho   = ncread(fname,'y_rho');
h       = ncread(fname,'h');

%% STEP 2: Create a z-grid by choosing N-levels
N       = 100;
[p,q]   = size(x_rho);

for i=1:1:p
    for j=1:1:q
        z(i,j,1:N)=-h(i,j):h(i,j)/(N-1):0;
    end
end

%% STEP 3: Choose a stratification profile and create temperature (Linear here)
dtdz = 0.25;
T0   = 20;
T1   = T0 + z*dtdz;

%% STEP 4: Get vertical stretching and z_grid in ROMS coordinates
Vtransform      = 2;
Vstretching     = 4;
theta_s         = 6;
theta_b         = 4;
Tcline          = 1;
N               = 10;
[s_w,Cs_w]      = stretching(Vstretching,theta_s,theta_b,Tcline,N,1);
[s_rho,Cs_rho]  = stretching(Vstretching,theta_s,theta_b,Tcline,N,0);
hc              = Tcline;
zeta            = 0*h;     %Assume no sea-surface elevation at t=0

[Hz,z_w,z_r]    = get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_rho,Cs_w,zeta,Vtransform);


%% STEP 5: Interpolate temperature to ROMS z-grid and define constant salinity
for i=1:1:p
    for j=1:1:q
        T2(i,j,:)= interp1(squeeze(z(i,j,:)),squeeze(T1(i,j,:)),squeeze(z_r(i,j,:)));
    end
end

temp     = T2;
salt     = 0*temp+35;

%% STEP 6: Plotting
figure
xx   = repmat(x_rho(:,1),[1 N]);
xx   = xx';
zz   = squeeze(z_r(:,1,:));
zz   = zz';
tt   = squeeze(temp(:,1,:));
tt   = tt';
pcolor(xx,zz,tt);
shading interp;
colormap('jet');
colorbar

%% STEP 7: Get rest of the variables:
zeta  = 0*h;
ubar  = 0*h;
vbar  = 0*h;
ubar  = ubar(1:end-1,:);
vbar  = vbar(:,end-1);
u     = 0*temp;
u     = u(1:end-1,:,:);
v     = 0*temp;
v     = v(:,1:end-1,:);
ocean_time = 0;

eval(['save ',sname,' temp',' salt',' u',' v',' ubar',' vbar',' zeta',' hc',' Tcline',' Vtransform',' Vstretching',' N',' theta_s', ' theta_b',' ocean_time']);

%% STEP 7: Create Initial Forcing File
create_roms_init_NK
