% --------------------------------------------------------------------------
% File Name: getval_savedata_mex.m
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
% getval_savedata_mex
% --------------------------------------------------------------------------


clear all
close all 
clc

% Set Seed
rng(15)

addpath savedata/

% Load parameters
P = parameters1;
P.bound = (P.sigma^2/(1-P.rho^2))^0.5;
O.delbound = [1-4*P.bound 1+4*P.bound];
O.del_pts = 21;
O.e_pts = 10;

% Chebyshev polynomial order in each dimension (must <= pts in grid)
O.n1 = 4;
value = zeros(length(P.pi_targ),1);

for i = 1:length(P.pi_targ)
    load(strcat('ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'_Ps',num2str(P.Ps),'_Pd',num2str(P.Pd),'.mat'))
    
    E_val = sim_val_zlb(P,O,C);
    value(i,1) = E_val;
    
    disp(strcat('Value calculated for PItarg =',num2str(400*(P.pi_targ(i)-1))))
end

max_value = max(value(:,1));
value_inx = find(max_value == value(:,1));
opt_inf = 400*(P.pi_targ(value_inx)-1);

opt_inf = [P.alpha,opt_inf];
disp(opt_inf);

fig(1) = figure(1);
box on
hold on
grid on
plot(400*(P.pi_targ-1),value,'k','LineWidth',2)
xlabel('\Pi^{targ}','FontSize',16)
title('Value','FontSize',16,'FontWeight','Normal')
set(gca,'XLim',[400*(P.pi_targ(1)-1) 400*(P.pi_targ(end)-1)],'FontSize',16)
line([opt_inf(2) opt_inf(2)],get(gca,'YLim'),'Color','b','LineStyle','--','LineWidth',1.5)


set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc',strcat('OptInfTarg','_Ps',num2str(P.Ps),'_Pd',num2str(P.Pd),'.eps'));