%--------------------------------------------------------------------------
% File Name: max_shock.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/iid/nstate/moments
% max_shock
%--------------------------------------------------------------------------

clear all
close all
clc

%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi_grid = [1.2, 2, 5];
cRstar = 1/400;
p_m = 0.5; %Probability of being in the middle state 
cSIGMAd_in_grid = zeros(length(cPHIpi_grid),1);

for k = 1:length(cPHIpi_grid)
    cPHIpi = cPHIpi_grid(k);
    cSHOCK_lb =((cRstar)*(cPHIpi - 1))/(cKAPPA*cPHIpi*cSIGMA);
    cSHOCK_ub =((cRstar)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi*cSIGMA);
    cSHOCK_bound = -(2*(cRstar)*(cPHIpi - 1)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi^2*cSIGMA*(cKAPPA*cSIGMA + 1)*(p_m - 1));
    
    cSIGMAd_in_grid(k) = max([cSHOCK_lb,cSHOCK_ub,cSHOCK_bound]);
end

nstate = 3;
% Find Middle State
rss_inx = (nstate+1)/2;

%% Main Code
for k = 1:length(cPHIpi_grid)
    cPHIpi = cPHIpi_grid(k);
    cSIGMAd_in = cSIGMAd_in_grid(k);
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd_in];
    
    options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'iter');

    x0 = cSIGMAd_in;
    func = @(x) param_solve_cSIGMAd(x,params,nstate,rss_inx);
    x_out_param = fsolve(func,x0,options);

    if k == 1;
        disp('Found maximum shock value consistent with equilibrium existence for low cPHPi case')
    elseif k == 2;
        disp('Found maximum shock value consistent with equilibrium existence for moderate cPHPi case')
    elseif k == 3;
        disp('Found maximum shock value consistent with equilibrium existence for high cPHPi case')
    end

    %% Save data
    if k == 1;
       cSIGMAd_max.l = x_out_param;
    elseif k == 2;
        cSIGMAd_max.m = x_out_param;
    elseif k == 3;
        cSIGMAd_max.h = x_out_param;
    end
    
end
save('cSIGMAd_max.mat','cSIGMAd_max')
