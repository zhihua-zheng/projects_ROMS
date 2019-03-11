%Clear Space
clear all %#ok<CLALL>

%% Define Coordinates
x       = -20:10:20;
y       = 0:10:40;
[X,Y]   = meshgrid(x,y);
X       = X'; % why transpose? row priority in ROMS, W-E as column.
Y       = Y';
h       = 0*X+4000;
mask_rho= h*0+1;
theta   = 50;
f       = 2*7.29*10^-5*sind(theta);
save uniform_depth_bathymetry_50m.mat X Y h mask_rho f