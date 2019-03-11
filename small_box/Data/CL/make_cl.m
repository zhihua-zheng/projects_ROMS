
%% processing climatology(cl)

load ts_1995_2017

z = sal_cl(1,:);
sal_cl = nanmean(sal_cl(2:5,:));
temp_cl = nanmean(temp_cl(2:5,:));
save('ocsp_cl_1995_2017.mat','z','sal_cl','temp_cl')
clear
