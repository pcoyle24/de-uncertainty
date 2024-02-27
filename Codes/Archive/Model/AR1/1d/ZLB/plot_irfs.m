% --------------------------------------------------------------------------
% File Name: plot_irfs.m
% Author: Philip Coyle
% Date Created: 11/09/2018
% 
% Add Appropriate Paths 
% codetoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/Code_Toolbox'; 
% mextoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/mex_functions';
% addpath(genpath(codetoolbox));
% addpath(genpath(mextoolbox)); 
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/AR1/1d/ZLB
% plot_irfs
% --------------------------------------------------------------------------


clear all
close all 
clc

% Set Seed
rng(22)

addpath savedata/
%% ------------------------------------------------------------------------
%  Initialize Policy Functions
%  ------------------------------------------------------------------------
% Load parameters
P = parameters;

%% ------------------------------------------------------------------------
%  IRF Analysis
%  ------------------------------------------------------------------------
time = 1:40;
per = 40;
c_irf = zeros(per,2);
pi_irf = zeros(per,2);
r_irf = zeros(per,2);

X = zeros(per,2,3);

shock = 4*P.bound;

for i = 1:length(P.pi_targ)
    load(strcat('ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'.mat'));
    
    c_init = S.c;
    pi_init = S.inf;
    r_init = S.r;
    del_yesterday = 1;
    
    for k = 1:per
        if k == 1
            del_today = P.rho*(del_yesterday - 1) + 1 + shock;
        else
            del_today = P.rho*(del_yesterday - 1) + 1;
        end
        
        r_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ar,C.max,C.T,C.P);
        if r_today >=1
            c_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac,C.max,C.T,C.P);
            pi_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf,C.max,C.T,C.P);
        else
            r_today = 1;
            c_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_zlb,C.max,C.T,C.P);
            pi_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_zlb,C.max,C.T,C.P);
        end
        
        c_irf(k,i) = 100*(c_today-1);
        pi_irf(k,i) = 400*(pi_today-1);
        r_irf(k,i) = 400*(r_today-1);
        
        del_yesterday = del_today;
    end
end 
X(:,:,1) = c_irf;
X(:,:,2) = pi_irf;
X(:,:,3) = r_irf;


%% ------------------------------------------------------------------------
%  Plot IRFs
%  ------------------------------------------------------------------------

title_fig = {'Consumption', 'Inflation Rate', 'Interest Rate'};
colors = {'k', 'b'};

fig(1) = figure(1);
for i = 1:length(title_fig)
    for k = 1:length(colors)
        subplot(2,3,i)
        box on
        hold on
        grid on
        if i == 1
            h(k) = plot(time,X(:,k,i),colors{k},'LineWidth',2);
        else
            plot(time,X(:,k,i),colors{k},'LineWidth',2)
        end
        xlabel('Periods','FontSize',16)
        title(title_fig{i},'FontSize',16,'FontWeight','Normal')
        set(gca,'XLim',[time(1) time(end)],'FontSize',16)
    end
end
L = legend([h(1) h(2)],'\Pi^{targ} = 1.5%', '\Pi^{targ} = 2%');
set(L,'Location','SouthEast','FontSize',15)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','IRFs.eps');