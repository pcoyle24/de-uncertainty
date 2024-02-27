% --------------------------------------------------------------------------
% File Name: plot_irfs_cPHIpi.m
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
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/AR1/1d/sunspot/RSSPFs
% plot_irfs_cPHIpi
% --------------------------------------------------------------------------



clear all
close all 
clc

addpath savedata/
%% ------------------------------------------------------------------------
%  Initialize Policy Functions
%  ------------------------------------------------------------------------
% Load parameters
P = parameters_cPHIpi;

%% ------------------------------------------------------------------------
%  IRF Analysis
%  ------------------------------------------------------------------------
time = 1:40;
per = 40;

c_d_irf = zeros(per,3);
pi_d_irf = zeros(per,3);
r_d_irf = zeros(per,3);

X = zeros(per,3,3);

for i = 1:length(P.phi_pi_grid)
    P.phi_pi = P.phi_pi_grid(i);
    shock = 4*P.bound;
       
    load(strcat('ChebPFs_cPHIpi',num2str(P.phi_pi),'.mat'));
    
    SS(1) = 100*(S.c_d - 1);
    SS(2) = 400*(S.inf_d-1);
    c_d_init = S.c_d;
    pi_d_init = S.inf_d;
    r_d_init = S.r_d;
    
    del_yesterday = 1;
    
    for k = 1:per
        if k == 1
            del_today = P.rho*(del_yesterday - 1) + 1 + shock;
        else
            del_today = P.rho*(del_yesterday - 1) + 1;
        end        
        % Deflationary Regime
        r_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ar_d,C.max,C.T,C.P);
        if r_d_today >=1
            c_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_d,C.max,C.T,C.P);
            pi_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_d,C.max,C.T,C.P);
        else
            r_d_today = 1;
            c_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_zlb_d,C.max,C.T,C.P);
            pi_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_zlb_d,C.max,C.T,C.P);
        end
        
        c_d_irf(k,i) = 100*(c_d_today-1);
        pi_d_irf(k,i) = 400*(pi_d_today-1);
        r_d_irf(k,i) = 400*(r_d_today-1);
        
        del_yesterday = del_today;
    end
end 
X(:,:,1) = c_d_irf;
X(:,:,2) = pi_d_irf;
X(:,:,3) = r_d_irf;


%% ------------------------------------------------------------------------
%  Plot IRFs
%  ------------------------------------------------------------------------

title_fig = {'Consumption (D)', 'Inflation Rate (D)', 'Interest Rate (D)'};
colors = {'k','b','r'};

fig(1) = figure(1);
for i = 1:length(title_fig)
    for k = 1:length(colors)
        subplot(2,3,i)
        box on
        hold on
        grid on
        if i == 3
            h(k) = plot(time,X(:,k,i),colors{k},'LineWidth',2);
        else
            plot(time,X(:,k,i),colors{k},'LineWidth',2)
        end
        if i == 1
            ylabel('%Dev from ESS','FontSize',16)
        else
            ylabel('% Annual Rate','FontSize',16)
        end
        if i == 1 || i == 2
            line(get(gca,'XLim'),[SS(i) SS(i)],'Color','k','LineStyle','-','LineWidth',0.5)
        end
        xlabel('Periods','FontSize',16)
        title(title_fig{i},'FontSize',16,'FontWeight','Normal')
        set(gca,'XLim',[time(1) time(end)],'FontSize',16)
    end
end
L = legend([h(1) h(2) h(3)],'\phi_{\pi} = 1.75','\phi_{\pi} = 3','\phi_{\pi} = 5');
set(L,'Location','SouthEast','FontSize',12)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
if P.alpha == 1
    print(fig(1),'-depsc','IRFs_SunSpot_VarycPHipi_cALPHA1.eps');
    print(fig(1),'-dpdf','IRFs_SunSpot_VarycPHipi_cALPHA1.pdf');
elseif P.alpha == 0
    print(fig(1),'-depsc','IRFs_SunSpot_VarycPHipi_cALPHA0.eps');
    print(fig(1),'-dpdf','IRFs_SunSpot_VarycPHipi_cALPHA0.pdf');
end