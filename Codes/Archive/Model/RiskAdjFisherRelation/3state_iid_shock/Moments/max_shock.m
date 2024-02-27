%--------------------------------------------------------------------------
% File Name: max_shock.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Codes/RiskAdjFisherRelation/3state_iid_shock/Moments
% max_shock
%--------------------------------------------------------------------------

clear all
close all
clc

%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 2;
cRstar = 1/400;
p_m = 0.5; %Probability of being in the middle state 

cSHOCK_lb =((cRstar)*(cPHIpi - 1))/(cKAPPA*cPHIpi*cSIGMA);
cSHOCK_ub =((cRstar)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi*cSIGMA);
cSHOCK_bound = -(2*(cRstar)*(cPHIpi - 1)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi^2*cSIGMA*(cKAPPA*cSIGMA + 1)*(p_m - 1));

cSIGMAd_in = max([cSHOCK_lb,cSHOCK_ub,cSHOCK_bound]);%- 0.000311840838465;
params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd_in];


nstate = 3;
% Find Middle State
rss_inx = (nstate+1)/2;

%% Main Code
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'iter');

x0 = cSIGMAd_in;
func = @(x) param_solve_cSIGMAd(x,params,nstate,rss_inx);
x_out_param = fsolve(func,x0,options);

disp('Found maximum shock value consistent with equilibrium existence')

%% Save data
params = x_out_param;
save('params.mat','params')
