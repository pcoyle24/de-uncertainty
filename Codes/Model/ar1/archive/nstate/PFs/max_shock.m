%--------------------------------------------------------------------------
% File Name: max_shock.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/ar1/nstate/PFs
% max_shock
%--------------------------------------------------------------------------

clear all
close all
clc

%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi_grid = [2, 4, 10];
cRstar = 1/400;
cRHO = 0.8; 

cSIGMAd_in_grid = [0.357/100, 0.367/100, 0.366/100];
for k = 1:length(cPHIpi_grid)
    cPHIpi = cPHIpi_grid(k);
    cSIGMAd_in = cSIGMAd_in_grid(k);
    
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd_in];
    nstate = 21;
    % Find Middle State
    rss_inx = (nstate+1)/2;

    %% Main Code
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
