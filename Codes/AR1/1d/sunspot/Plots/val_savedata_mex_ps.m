% --------------------------------------------------------------------------
% File Name: getval_savedata_mex_ps_pd.m
% Author: Philip Coyle
% Date Created: 02/06/2018
% Last Updated: 07/19/2018
% 
% Add Appropriate Paths %
% codetoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/Code_Toolbox'; 
% mextoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/mex_functions';
% addpath(genpath(codetoolbox));
% addpath(genpath(mextoolbox)); 
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/AR1/1d/sunspot
% getval_savedata_mex_ps_pd
% --------------------------------------------------------------------------


clear all
close all 
clc

% Set Seed
% rng(22)
% rng(42)
rng(15)

addpath savedata/

% Load parameters
P = parameters2;
P.bound = (P.sigma^2/(1-P.rho^2))^0.5;
O.delbound = [1-4*P.bound 1+4*P.bound];
O.del_pts = 21;
O.e_pts = 10;

% Chebyshev polynomial order in each dimension (must <= pts in grid)
O.n1 = 4;
value = zeros(length(P.pi_targ),length(P.Ps_grid));
for s = length(P.Ps_grid):-1:1
    P.Ps = P.Ps_grid(s);
    disp(strcat('P_s =',num2str(P.Ps))); 
    for i = 1:length(P.pi_targ)
        load(strcat('ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'_Ps',num2str(P.Ps),'_Pd',num2str(P.Pd),'.mat'))

        E_val = sim_val_zlb(P,O,C);
        value(i,s) = E_val;

        disp(strcat('Value calculated for PItarg =',num2str(400*(P.pi_targ(i)-1))))
    end
end

opt_inf = zeros(length(P.Ps_grid),1);
for s = 1:length(P.Ps_grid)
    max_value = max(value(:,s));
    value_inx = find(max_value == value(:,s));
    opt_inf(s,1) = 400*(P.pi_targ(value_inx)-1);
end

opt_inf_disp = [P.Ps_grid',opt_inf];
disp(opt_inf_disp);

fig(1) = figure(1);
box on
hold on
grid on
plot(P.Ps_grid,opt_inf,'k','LineWidth',2)
xlabel('p_T','FontSize',16)
xlabel('Optimal Inflation Target','FontSize',16)
% title('O','FontSize',16,'FontWeight','Normal')
set(gca,'XLim',[P.Ps_grid(1) P.Ps_grid(end)],'FontSize',16)


set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc',strcat('OptInfTarg_VaryPsPd.eps'));