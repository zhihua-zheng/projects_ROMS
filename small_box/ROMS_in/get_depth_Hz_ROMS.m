function [Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta,Vtransform)
% function [Hz,z_w,z_r]=get_depth_Hz_ROMS(h,hc,s_rho,s_w,Cs_r,Cs_w,zeta,Vtransform);
% This function computes the vertical grid metrics, z_r,z_w and Hz for a
% ROMS grid as a function of time:
% INPUT: 
% h    = Water depth (+ve in meters)
% hc   = Critical depth controlling the vertical stretching (in m)
% s_rho= S coordinates at RHO points
% s_w  = S coordinates at w points
% Cs_r = S coordinate stretching curves at rho points
% Cs_w = S coordinate stretching curves at w points
% zeta = sea surface elevation (m)
% Vtransform = 1 should be default?
% OUTPUT
% Hz(1:Lp,1:Mp,1:N), z_r(1:Lp,1:Mp,1:N), z_w(1:Lp,1:Mp,1:N) 
% 

N1         = length(s_rho);
N2         = length(s_w);
j          = h==0;
h(j)       = eps;
z_w(1,:,:) = -h;

cff_r      = hc*s_rho;
cff_w      = hc*s_w;
cff1_r     = Cs_r;
cff1_w     = Cs_w;
hinv       = 1./(hc+h);
hwater     = h;
count=0;

for i=1:1:N1
    count = count+1;
    if Vtransform==1
        cff2_r= (cff_r(i)-cff1_r(i)*hc)+h*cff1_r(i);
        cff2_w= (cff_w(i+1)-cff1_w(i+1)*hc)+h*cff1_w(i+1);
        z_w(i+1,:,:) = cff2_w + zeta.*(1+cff2_w./hwater);
        z_r(i,:,:)   = cff2_r + zeta.*(1+cff2_r./hwater);
        Hz(i,:,:)    = z_w(i+1,:,:)-z_w(i,:,:);
    elseif Vtransform==2
        cff2_r= (cff_r(i)+cff1_r(i)*hwater).*hinv;
        cff2_w= (cff_w(i+1)+cff1_w(i+1)*hwater).*hinv;
        z_w(i+1,:,:) = zeta + (zeta+hwater).*cff2_w;
        z_r(i,:,:)   = zeta + (zeta+hwater).*cff2_r;
        Hz(i,:,:)    = z_w(i+1,:,:)-z_w(i,:,:);
    end
        
end

Hz    = permute(Hz,[2 3 1]);
z_r   = permute(z_r,[2 3 1]);
z_w   = permute(z_w,[2 3 1]);
