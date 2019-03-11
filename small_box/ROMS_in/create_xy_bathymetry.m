clear all %#ok<CLALL>
% addpath(genpath('/home/nirni/Models/COAWST/Tools/mfiles'));

%% STEP 1: Load Bathymtery File
load uniform_depth_bathymetry_50m.mat
spherical   = 'F';
projection  ='mercator';

%% STEP 2: Provide grid name
fname        ='uniform_grid_50m.nc';
%%%%%%%%%%%%% END OF USER SECTION %%%%%%%%%%%%%%%%

%create roms grid
  rho.x    = X;
  rho.y    = Y;  
  rho.depth= h;
  rho.mask = mask_rho;
  rho.f    = f;
  save temp_jcw33.mat
  
  % double the apostrophe 
  eval(['mat2roms_mw_nk(''temp_jcw33.mat'',''', fname, ''');'])
  !rm temp_jcw33.mat
  disp(['Created roms grid -->   ',fname])

%to create swan grid then use:
  disp('Created swan grid + bathy files:')
  roms2swan(rho.x,rho.y,rho.depth,rho.mask);


