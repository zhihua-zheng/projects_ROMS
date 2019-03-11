%This function creates an initial forcing file for the Regional Ocean Modeling System 
% using an input of prognostic variables u,v,salt,temp,ubar,vbar and zeta
% from the parent grid.
clear 
% addpath(genpath('/home/n2kumar/MAT_Mytoolbox/MATLAB_ROMS/Utility'));
%% Step 0: Input stretching and vertical co-ordinate parameters
Vstretching     = 2;
theta_s         = 6;
theta_b         = 2.5;
Tcline          = 1; %Nirni changed this from 10 to 1 to -2
N               = 15;
[s_w,Cs_w]      = stretching(Vstretching,theta_s,theta_b,Tcline,N,1,0);
[s_rho,Cs_rho]  = stretching(Vstretching,theta_s,theta_b,Tcline,N,0,0);


%% Step 1:Load the parent grid interpolated initial contiions:
%A               = load('Init_Forcing_MAT/inital_forcing_600_grid_0006.mat');
%A               = load('inital_forcing_22_grid_Jun_nw.mat');
u               = A.u;
v               = A.v;
salt            = A.salt;
temp            = A.temp;
ubar            = A.ubar;
vbar            = A.vbar;
zeta            = A.zeta;

%% Step 2: Load child grid and init file name
cgrid           = '../Grids/GT_Grids/interp-22m-grid_v3.nc';
lon_rho         =  ncread(cgrid,'lon_rho');
h               =  ncread(cgrid,'h');
fn              = 'init_forcing_interp-22m-grid_jun_nw.nc';

%% Step 3: Create time stamp
tstart1      = datenum(1999,1,1,0,0,0);
tstart       = datenum(2000,5,30,12,0,0);
ocean_time   = tstart-tstart1;

%% Step 4: Parse Data
[ix,iy]     = size(zeta);
[LP,MP]     = size(lon_rho);
L       = LP-1;
Lm      = L-1;
M       = MP-1;
Mm      = M-1;
xpsi    = L;
xrho    = LP;
xu      = L;
xv      = LP;
epsi    = M;
erho    = MP;
eu      = MP;
ev      = M;
s       = N;
hc      = min(min(h));

if ~isequal(hc,Tcline)
    hc=Tcline;
end

%% Step 6:
NAT     = 2;  % Number of active tracers. Usually 2 (temp + salt). Same 

%% Step 7: Enter the number of mud sediments (NCS) and number of sand sediments (NNS).
NCS     = 0;   %number of cohesive sed classes
NNS     = 1;   %number of non-cohesive sed classes

% Calculate sediment parameters. Do not alter.
NST     = NCS + NNS;     % total number of sed tracers.
NT      = NAT+NST;       % total number of tracers.

% Enter number of bed layers (This value should be the same as in
%mod_param.F)
Nbed    = 1;

% Sediment class properties (in order, mud first then sand).
% These values should coincide with your sediment.in file.
mud_Srho    = ones(1,NCS)*2650;     %kg m-3, NCS values
mud_Sd50    = [0.01 0.06]/1000;     %m,      NCS values
mud_Wsed    = [0.10 1.0]/1000;      %m s-1,  NCS values
mud_tau_ce  = [0.05 0.05];          %N m-2,  NCS values
mud_Erate   = [5 5]*1e-5;           %kg m-2 s-1, NCS values
sand_Srho   = ones(1,NNS)*2650;     %kg m-3, NNS values
sand_Sd50   = [0.08]/1000;           %m,      NNS values
sand_Wsed   = [1.0]/1000;           %m s-1,  NNS values
sand_tau_ce = [0.07];               %N m-2,  NNS values
sand_Erate  = [1]*1e-4;             %kg m-2 s-1, NNS values

Srho        = [mud_Srho,sand_Srho];
Sd50        = [mud_Sd50,sand_Sd50];
Wsed        = [mud_Wsed,sand_Wsed];
tau_ce      = [mud_tau_ce,sand_tau_ce];
Erate       = [mud_Erate,sand_Erate];

% Provide initial sediment properties in water column.

display('Initializing suspended sediments.')
%
% mud.
%
for idmud=1:NCS
  count=['0',num2str(idmud)];
  count=count(end-1:end);
  eval(['mud_',count,'(1:xi_rho,1:eta_rho,1:N,1:length(ocean_time)) = 0;'])   %mud conc in water column
end
%
% sand.
%
for isand=1:NNS
  count=['0',num2str(isand)];
  count=count(end-1:end);
  eval(['sand_',count,'(1:xrho,1:erho,1:N,1:length(ocean_time)) = 0;'])  %sand conc in water column
end

% Provide initial sediment properties in bed.
%
% bed properties
  display('Initializing sediment bed.')
%
for k=1:Nbed
  for time=1:length(ocean_time)
    bed_thickness(1:xrho,1:erho,k,time) = 0.1;
    bed_age(1:xrho,1:erho,k,time)       = ocean_time(1);
    bed_porosity(1:xrho,1:erho,k,time)  = 0.9;
    bed_biodiff(1:xrho,1:erho,k,time)   = 0.0;
  end
end
%
% for mud
%
for idsed=1:NCS
  count=['0',num2str(idsed)];
  count=count(end-1:end);
  for k=1:Nbed
    for time=1:length(ocean_time)
      eval(['mudfrac_',count,'(1:xrho,1:erho,k,time) = 1/NST;'])      %fraction of each sed class in each bed cell
      eval(['mudmass_',count,'(1:xrho,1:erho,k,time) = squeeze(bed_thickness(1:xrho,1:erho,k,time)).*Srho(idsed).*(1.0-squeeze(bed_porosity(1:xrho,1:erho,k,time))).*squeeze(mudfrac_',count,'(1:xrho,1:erho,k,time));'])          %mass of each sed class in each bed cell
    end
  end
end
%
% for sand
%
for isand=1:NNS
 count=['0',num2str(isand)];
 count=count(end-1:end);
  for k=1:Nbed
    for time=1:length(ocean_time)
      eval(['sandfrac_',count,'(1:xrho,1:erho,k,time) = 1/NST;'])      %fraction of each sed class in each bed cell
      eval(['sandmass_',count,'(1:xrho,1:erho,k,time) = squeeze(bed_thickness(1:xrho,1:erho,k,time)).*Srho(isand).*(1.0-squeeze(bed_porosity(1:xrho,1:erho,k,time))).*squeeze(sandfrac_',count,'(1:xrho,1:erho,k,time));'])          %mass of each sed class in each bed cell
    end
  end
end

%
% set some surface properties
%
display('Initializing sediment surface properties.')
%
cff1=1.0;
cff2=1.0;
cff3=1.0;
cff4=1.0;
for ised=1:NCS
    count=['0',num2str(ised)];
    count=count(end-1:end);
    eval(['cff1=cff1*mud_Sd50(ised)^squeeze(mudfrac_',count,'(1:xrho,1:erho,1,1));'])
    eval(['cff2=cff2*mud_Srho(ised)^squeeze(mudfrac_',count,'(1:xrho,1:erho,1,1));'])
    eval(['cff3=cff3*mud_Wsed(ised)^squeeze(mudfrac_',count,'(1:xrho,1:erho,1,1));'])
    eval(['cff4=cff4*mud_tau_ce(ised)^squeeze(mudfrac_',count,'(1:xrho,1:erho,1,1));'])
end
for ised=1:NNS
    count=['0',num2str(ised)];
    count=count(end-1:end);
    eval(['cff1=cff1*sand_Sd50(ised).^squeeze(sandfrac_',count,'(1:xrho,1:erho,1,1));'])
    eval(['cff2=cff2*sand_Srho(ised).^squeeze(sandfrac_',count,'(1:xrho,1:erho,1,1));'])
    eval(['cff3=cff3*sand_Wsed(ised).^squeeze(sandfrac_',count,'(1:xrho,1:erho,1,1));'])
    eval(['cff4=cff4*sand_tau_ce(ised).^squeeze(sandfrac_',count,'(1:xrho,1:erho,1,1));'])
end
grain_diameter(1:xrho,1:erho,time) = cff1;
grain_density(1:xrho,1:erho,time)  = cff2;
settling_vel(1:xrho,1:erho,time)   = cff3;
erosion_stress(1:xrho,1:erho,time) = cff4;
ripple_length(1:xrho,1:erho,time)  = 0.10;
ripple_height(1:xrho,1:erho,time)  = 0.01;
dmix_offset(1:xrho,1:erho,time)    = 0.0;
dmix_slope(1:xrho,1:erho,time)     = 0.0;
dmix_time(1:xrho,1:erho,time)      = 0.0;

nc_init=netcdf.create(fn,'clobber');
if isempty(nc_init), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by' mfilename 'on ' datestr(now)]);
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'type', 'initial forcing file from a bigger grid');

%% Dimensions:
disp(' ## Defining Dimensions...')
 
psidimID   = netcdf.defDim(nc_init,'xi_psi',L);
xrhodimID  = netcdf.defDim(nc_init,'xi_rho',LP);
xudimID    = netcdf.defDim(nc_init,'xi_u',L);
xvdimID    = netcdf.defDim(nc_init,'xi_v',LP);

epsidimID  = netcdf.defDim(nc_init,'eta_psi',M);
erhodimID  = netcdf.defDim(nc_init,'eta_rho',MP);
eudimID    = netcdf.defDim(nc_init,'eta_u',MP);
evdimID    = netcdf.defDim(nc_init,'eta_v',M);
s_rhodimID = netcdf.defDim(nc_init,'s_rho',length(s_rho));
s_wdimID   = netcdf.defDim(nc_init,'s_w',length(s_w));
timedimID  = netcdf.defDim(nc_init,'temp_time',length(ocean_time));
oneID      = netcdf.defDim(nc_init,'one',1);
NbeddimID  = netcdf.defDim(nc_init,'Nbed',Nbed);
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

sphericalID = netcdf.defVar(nc_init,'spherical','char',oneID);
netcdf.putAtt(nc_init,sphericalID,'long_name','grid type logical switch');
netcdf.putAtt(nc_init,sphericalID,'flag_meanings','spherical Cartesian');
netcdf.putAtt(nc_init,sphericalID,'flag_values','T, F');

VtransformID = netcdf.defVar(nc_init,'Vtransform','long',oneID);
netcdf.putAtt(nc_init,VtransformID,'long_name','vertical terrain-following transformation equation');

VstretchingID = netcdf.defVar(nc_init,'Vstretching','long',oneID);
netcdf.putAtt(nc_init,VstretchingID,'long_name','vertical terrain-following stretching function');
 
thetasID = netcdf.defVar(nc_init,'theta_s','double',oneID);
netcdf.putAtt(nc_init,thetasID,'long_name','S-coordinate surface control parameter');
netcdf.putAtt(nc_init,thetasID,'units','nondimensional');
netcdf.putAtt(nc_init,thetasID,'field','theta_s, scalar, series');

thetabID = netcdf.defVar(nc_init,'theta_b','double',oneID);
netcdf.putAtt(nc_init,thetabID,'long_name','S-coordinate bottom control parameter');
netcdf.putAtt(nc_init,thetabID,'units','nondimensional');
netcdf.putAtt(nc_init,thetabID,'field','theta_b, scalar, series');

tclineID = netcdf.defVar(nc_init,'Tcline','double',oneID);
netcdf.putAtt(nc_init,tclineID,'long_name','S-coordinate surface/bottom layer width');
netcdf.putAtt(nc_init,tclineID,'units','meter');
netcdf.putAtt(nc_init,tclineID,'field','Tcline, scalar, series');

hcID = netcdf.defVar(nc_init,'hc','double',oneID);
netcdf.putAtt(nc_init,hcID,'long_name','S-coordinate parameter, critical depth');
netcdf.putAtt(nc_init,hcID,'units','meter');
netcdf.putAtt(nc_init,hcID,'field','hc, scalar, series');

csrID = netcdf.defVar(nc_init,'Cs_r','double',s_rhodimID);
netcdf.putAtt(nc_init,csrID,'long_name','S-coordinate stretching curves at RHO-points');
netcdf.putAtt(nc_init,csrID,'units','nondimensional');
netcdf.putAtt(nc_init,csrID,'valid_min',-1);
netcdf.putAtt(nc_init,csrID,'valid_max',0);
netcdf.putAtt(nc_init,csrID,'field','Cs_r, scalar, series');

scrID = netcdf.defVar(nc_init,'sc_r','double',s_rhodimID);
netcdf.putAtt(nc_init,scrID,'long_name','S-coordinate at RHO-points');
netcdf.putAtt(nc_init,scrID,'units','nondimensional');
netcdf.putAtt(nc_init,scrID,'valid_min',-1);
netcdf.putAtt(nc_init,scrID,'valid_max',0);
netcdf.putAtt(nc_init,scrID,'field','sc_r, scalar, series');

cswID = netcdf.defVar(nc_init,'Cs_w','double',s_wdimID);
netcdf.putAtt(nc_init,cswID,'long_name','S-coordinate stretching curves at W-points');
netcdf.putAtt(nc_init,cswID,'units','nondimensional');
netcdf.putAtt(nc_init,cswID,'valid_min',-1);
netcdf.putAtt(nc_init,cswID,'valid_max',0);
netcdf.putAtt(nc_init,cswID,'field','Cs_w, scalar, series');

scwID = netcdf.defVar(nc_init,'sc_w','double',s_wdimID);
netcdf.putAtt(nc_init,scwID,'long_name','S-coordinate at W-points');
netcdf.putAtt(nc_init,scwID,'units','nondimensional');
netcdf.putAtt(nc_init,scwID,'field','sc_w, scalar, series');
netcdf.putAtt(nc_init,scwID,'valid_min',-1);
netcdf.putAtt(nc_init,scwID,'valid_max',0);

ocean_timeID = netcdf.defVar(nc_init,'ocean_time','double',timedimID);
netcdf.putAtt(nc_init,ocean_timeID,'long_name','time since initialization');
netcdf.putAtt(nc_init,ocean_timeID,'units','days');
netcdf.putAtt(nc_init,ocean_timeID,'field','ocean_time, scalar, series');

saltID = netcdf.defVar(nc_init,'salt','float',[xrhodimID erhodimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,saltID,'long_name','salinity');
netcdf.putAtt(nc_init,saltID,'units','PSU');
netcdf.putAtt(nc_init,saltID,'field','salinity, scalar, series');

tempID = netcdf.defVar(nc_init,'temp','float',[xrhodimID erhodimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,tempID,'long_name','temperature');
netcdf.putAtt(nc_init,tempID,'units','C');
netcdf.putAtt(nc_init,tempID,'field','temperature, scalar, series');

uID = netcdf.defVar(nc_init,'u','float',[xudimID eudimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,uID,'long_name','u-momentum component');
netcdf.putAtt(nc_init,uID,'units','meter second-1');
netcdf.putAtt(nc_init,uID,'field','u-velocity, scalar, series');

ubarID = netcdf.defVar(nc_init,'ubar','float',[xudimID eudimID timedimID]);
netcdf.putAtt(nc_init,ubarID,'long_name','vertically integrated u-momentum component');
netcdf.putAtt(nc_init,ubarID,'units','meter second-1');
netcdf.putAtt(nc_init,ubarID,'field','ubar-velocity, scalar, series');

vID = netcdf.defVar(nc_init,'v','float',[xvdimID evdimID s_rhodimID timedimID]);
netcdf.putAtt(nc_init,vID,'long_name','v-momentum component');
netcdf.putAtt(nc_init,vID,'units','meter second-1');
netcdf.putAtt(nc_init,vID,'field','v-velocity, scalar, series');

vbarID = netcdf.defVar(nc_init,'vbar','float',[xvdimID evdimID timedimID]);
netcdf.putAtt(nc_init,vbarID,'long_name','vertically integrated v-momentum component');
netcdf.putAtt(nc_init,vbarID,'units','meter second-1');
netcdf.putAtt(nc_init,vbarID,'field','vbar-velocity, scalar, series');
 
zetaID = netcdf.defVar(nc_init,'zeta','float',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,zetaID,'long_name','free-surface');
netcdf.putAtt(nc_init,zetaID,'units','meter');
netcdf.putAtt(nc_init,zetaID,'field','free-surface, scalar, series');

for mm=1:NCS
    count=['00',num2str(mm)];
    count=count(end-1:end);

    eval(['mud_',count,'ID = netcdf.defVar(nc_init,''mud_',count,''',''double'',[xrhodimID erhodimID s_rhodimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''long_name'',''suspended cohesive sediment, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''units'',''kilogram meter-3'');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mud_',count,'ID,''field'',''mud_',count,', scalar, series'');'])

    eval(['mudfrac_',count,'ID = netcdf.defVar(nc_init,''mudfrac_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''long_name'',''cohesive sediment fraction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''units'',''nondimensional'');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mudfrac_',count,'ID,''field'',''mudfrac_',count,', scalar, series'');'])

    eval(['mudmass_',count,'ID = netcdf.defVar(nc_init,''mudmass_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''long_name'',''cohesive sediment mass, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''units'',''kilogram meter-2'');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,mudmass_',count,'ID,''field'',''mudmass_',count,', scalar, series'');'])

end
for mm=1:NNS
    count=['00',num2str(mm)];
    count=count(end-1:end);

    eval(['sand_',count,'ID = netcdf.defVar(nc_init,''sand_',count,''',''double'',[xrhodimID erhodimID s_rhodimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''long_name'',''suspended noncohesive sediment, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''units'',''kilogram meter-3'');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sand_',count,'ID,''field'',''sand_',count,', scalar, series'');'])

    eval(['sandfrac_',count,'ID = netcdf.defVar(nc_init,''sandfrac_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''long_name'',''noncohesive sediment fraction, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''units'',''nondimensional'');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sandfrac_',count,'ID,''field'',''sandfrac_',count,', scalar, series'');'])

    eval(['sandmass_',count,'ID = netcdf.defVar(nc_init,''sandmass_',count,''',''double'',[xrhodimID erhodimID NbeddimID timedimID]);'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''long_name'',''noncohesive sediment mass, size class ',count,''');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''units'',''kilogram meter-2'');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''time'',''ocean_time'');'])
    eval(['netcdf.putAtt(nc_init,sandmass_',count,'ID,''field'',''sandmass_',count,', scalar, series'');'])

end

bed_thicknessID = netcdf.defVar(nc_init,'bed_thickness','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_thicknessID,'long_name','sediment layer thickness');
netcdf.putAtt(nc_init,bed_thicknessID,'units','meter');
netcdf.putAtt(nc_init,bed_thicknessID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_thicknessID,'field','bed thickness, scalar, series');

bed_ageID = netcdf.defVar(nc_init,'bed_age','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_ageID,'long_name','sediment layer age');
netcdf.putAtt(nc_init,bed_ageID,'units','day');
netcdf.putAtt(nc_init,bed_ageID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_ageID,'field','bed age, scalar, series');

bed_porosityID = netcdf.defVar(nc_init,'bed_porosity','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_porosityID,'long_name','sediment layer porosity');
netcdf.putAtt(nc_init,bed_porosityID,'units','nondimensional');
netcdf.putAtt(nc_init,bed_porosityID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_porosityID,'field','bed porosity, scalar, series');

bed_biodiffID = netcdf.defVar(nc_init,'bed_biodiff','double',[xrhodimID erhodimID NbeddimID timedimID]);
netcdf.putAtt(nc_init,bed_biodiffID,'long_name','biodiffusivity at bottom of each layer');
netcdf.putAtt(nc_init,bed_biodiffID,'units','meter2 second-1');
netcdf.putAtt(nc_init,bed_biodiffID,'time','ocean_time');
netcdf.putAtt(nc_init,bed_biodiffID,'field','bed biodiffusivity, scalar, series');

grain_diameterID = netcdf.defVar(nc_init,'grain_diameter','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,grain_diameterID,'long_name','sediment median grain diameter size');
netcdf.putAtt(nc_init,grain_diameterID,'units','meter');
netcdf.putAtt(nc_init,grain_diameterID,'time','ocean_time');
netcdf.putAtt(nc_init,grain_diameterID,'field','grain diameter, scalar, series');

grain_densityID = netcdf.defVar(nc_init,'grain_density','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,grain_densityID,'long_name','sediment median grain density');
netcdf.putAtt(nc_init,grain_densityID,'units','kilogram meter-3');
netcdf.putAtt(nc_init,grain_densityID,'time','ocean_time');
netcdf.putAtt(nc_init,grain_densityID,'field','grain density, scalar, series');

settling_velID = netcdf.defVar(nc_init,'settling_vel','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,settling_velID,'long_name','sediment median grain settling velocity');
netcdf.putAtt(nc_init,settling_velID,'units','meter second-1');
netcdf.putAtt(nc_init,settling_velID,'time','ocean_time');
netcdf.putAtt(nc_init,settling_velID,'field','settling vel, scalar, series');

erosion_stressID = netcdf.defVar(nc_init,'erosion_stress','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,erosion_stressID,'long_name','sediment median critical erosion stress');
netcdf.putAtt(nc_init,erosion_stressID,'units','meter2 second-2');
netcdf.putAtt(nc_init,erosion_stressID,'time','ocean_time');
netcdf.putAtt(nc_init,erosion_stressID,'field','erosion stress, scalar, series');

ripple_lengthID = netcdf.defVar(nc_init,'ripple_length','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,ripple_lengthID,'long_name','bottom ripple length');
netcdf.putAtt(nc_init,ripple_lengthID,'units','meter');
netcdf.putAtt(nc_init,ripple_lengthID,'time','ocean_time');
netcdf.putAtt(nc_init,ripple_lengthID,'field','ripple length, scalar, series');

ripple_heightID = netcdf.defVar(nc_init,'ripple_height','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,ripple_heightID,'long_name','bottom ripple height');
netcdf.putAtt(nc_init,ripple_heightID,'units','meter');
netcdf.putAtt(nc_init,ripple_heightID,'time','ocean_time');
netcdf.putAtt(nc_init,ripple_heightID,'field','ripple height, scalar, series');

dmix_offsetID = netcdf.defVar(nc_init,'dmix_offset','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,dmix_offsetID,'long_name','dmix erodibility profile offset');
netcdf.putAtt(nc_init,dmix_offsetID,'units','meter');
netcdf.putAtt(nc_init,dmix_offsetID,'time','ocean_time');
netcdf.putAtt(nc_init,dmix_offsetID,'field','dmix_offset, scalar, series');

dmix_slopeID = netcdf.defVar(nc_init,'dmix_slope','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,dmix_slopeID,'long_name','dmix erodibility profile slope');
netcdf.putAtt(nc_init,dmix_slopeID,'units','_');
netcdf.putAtt(nc_init,dmix_slopeID,'time','ocean_time');
netcdf.putAtt(nc_init,dmix_slopeID,'field','dmix_slope, scalar, series');

dmix_timeID = netcdf.defVar(nc_init,'dmix_time','double',[xrhodimID erhodimID timedimID]);
netcdf.putAtt(nc_init,dmix_timeID,'long_name','dmix erodibility profile time scale');
netcdf.putAtt(nc_init,dmix_timeID,'units','seconds');
netcdf.putAtt(nc_init,dmix_timeID,'time','ocean_time');
netcdf.putAtt(nc_init,dmix_timeID,'field','dmix_time, scalar, series');

%close file
netcdf.close(nc_init)

nc_init=netcdf.open(fn,'NC_WRITE');
%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

%% Vertical Co-ordinate Parameters

tempid = netcdf.inqVarID(nc_init,'Vtransform');
netcdf.putVar(nc_init,tempid,1);

tempid = netcdf.inqVarID(nc_init,'Vstretching');
netcdf.putVar(nc_init,tempid,1);

sphericalfid = netcdf.inqVarID(nc_init,'spherical');
netcdf.putVar(nc_init,sphericalfid,'T');

thetasid = netcdf.inqVarID(nc_init,'theta_s');
netcdf.putVar(nc_init,thetasid,theta_s);

thetabid = netcdf.inqVarID(nc_init,'theta_b');
netcdf.putVar(nc_init,thetabid,theta_b);

tclineid = netcdf.inqVarID(nc_init,'Tcline');
netcdf.putVar(nc_init,tclineid,Tcline);

hcid     = netcdf.inqVarID(nc_init,'hc');
netcdf.putVar(nc_init,hcid,hc);

csrid    = netcdf.inqVarID(nc_init,'Cs_r');
netcdf.putVar(nc_init,csrid,Cs_rho);

cswid    = netcdf.inqVarID(nc_init,'Cs_w');
netcdf.putVar(nc_init,cswid,Cs_w);

scrid    = netcdf.inqVarID(nc_init,'sc_r');
netcdf.putVar(nc_init,scrid,s_rho);

scwrid   = netcdf.inqVarID(nc_init,'sc_w');
netcdf.putVar(nc_init,scwrid,s_w);

%% Time
tempid = netcdf.inqVarID(nc_init,'ocean_time');  %get id
netcdf.putVar(nc_init,tempid,ocean_time);        %set variable

%% zeta , ubar, vbar
vars2d={'zeta','ubar','vbar'};
for i=1:length(vars2d)
    eval(['tempid = netcdf.inqVarID(nc_init,''',vars2d{i},''');']);%get id
    eval(['netcdf.putVar(nc_init,tempid,',vars2d{i},');']);%set variable
end
%% u,v, temp and salt
vars3d={'u','v','temp','salt'};
for i=1:length(vars3d)
    eval(['tempid = netcdf.inqVarID(nc_init,''',vars3d{i},''');']);%get id
    eval(['netcdf.putVar(nc_init,tempid,',vars3d{i},');']);%set variable
end

%% sediment parameters
tempid = netcdf.inqVarID(nc_init,'bed_thickness');  %get id
netcdf.putVar(nc_init,tempid,bed_thickness);        %set variable

tempid = netcdf.inqVarID(nc_init,'bed_age');        %get id
netcdf.putVar(nc_init,tempid,bed_thickness);        %set variable

tempid = netcdf.inqVarID(nc_init,'bed_porosity');   %get id
netcdf.putVar(nc_init,tempid,bed_thickness);        %set variable

tempid = netcdf.inqVarID(nc_init,'bed_biodiff');    %get id
netcdf.putVar(nc_init,tempid,bed_biodiff);          %set variable

tempid = netcdf.inqVarID(nc_init,'grain_diameter');  %get id
netcdf.putVar(nc_init,tempid,grain_diameter);        %set variable

tempid = netcdf.inqVarID(nc_init,'grain_density');  %get id
netcdf.putVar(nc_init,tempid,grain_density);        %set variable

tempid = netcdf.inqVarID(nc_init,'settling_vel');  %get id
netcdf.putVar(nc_init,tempid,settling_vel);        %set variable

tempid = netcdf.inqVarID(nc_init,'erosion_stress');  %get id
netcdf.putVar(nc_init,tempid,erosion_stress);        %set variable

tempid = netcdf.inqVarID(nc_init,'ripple_length');  %get id
netcdf.putVar(nc_init,tempid,ripple_length);        %set variable

tempid = netcdf.inqVarID(nc_init,'ripple_height');  %get id
netcdf.putVar(nc_init,tempid,ripple_height);        %set variable

tempid = netcdf.inqVarID(nc_init,'dmix_offset');  %get id
netcdf.putVar(nc_init,tempid,dmix_offset);        %set variable

tempid = netcdf.inqVarID(nc_init,'dmix_slope');  %get id
netcdf.putVar(nc_init,tempid,dmix_slope);        %set variable

tempid = netcdf.inqVarID(nc_init,'dmix_time');  %get id
netcdf.putVar(nc_init,tempid,dmix_time);        %set variable

for mm=1:NCS
    count=['00',num2str(mm)];
    count=count(end-1:end);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''mud_',count,''');']);%sed conc in water column
    eval(['netcdf.putVar(nc_init,tempid,mud_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''mudfrac_',count,''');']); %fraction of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,mudfrac_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''mudmass_',count,''');']);%mass of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,mudmass_',count,');']);
end
for mm=1:NNS
    count=['00',num2str(mm)];
    count=count(end-1:end);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''sand_',count,''');']);%sed conc in water column
    eval(['netcdf.putVar(nc_init,tempid,sand_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''sandfrac_',count,''');']); %fraction of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,sandfrac_',count,');']);
    
    eval(['tempid = netcdf.inqVarID(nc_init,''sandmass_',count,''');']);%mass of each sed class in each bed cell
    eval(['netcdf.putVar(nc_init,tempid,sandmass_',count,');']);
end


%%
netcdf.close(nc_init);
