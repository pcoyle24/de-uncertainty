%--------------------------------------------------------------------------
% File Name: stylized_iid_moments_cPHIpi_l.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/iid/nstate/moments
% stylized_iid_moments_cPHIpi_l
%--------------------------------------------------------------------------

clear all
close all
clc

load ('cSIGMAd_max.mat')

%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 1.2;
cRstar = 1/400;
p_m = 0.5; %Probability of being in the middle state 

cSIGMAd_grid = linspace(0,cSIGMAd_max.l,51);
% cSIGMAd_grid = cSIGMAd_max;

nstate = 3;
% Find Middle State
rss_inx = (nstate+1)/2;

%% Housekeeping
e_y_tr = zeros(length(cSIGMAd_grid),1);
e_y_dr = zeros(length(cSIGMAd_grid),1);

rss_y_tr = zeros(length(cSIGMAd_grid),1);
rss_y_dr = zeros(length(cSIGMAd_grid),1);

e_pi_tr = zeros(length(cSIGMAd_grid),1);
e_pi_dr = zeros(length(cSIGMAd_grid),1);

rss_pi_tr = zeros(length(cSIGMAd_grid),1);
rss_pi_dr = zeros(length(cSIGMAd_grid),1);

e_i_tr = zeros(length(cSIGMAd_grid),1);
e_i_dr = zeros(length(cSIGMAd_grid),1);

rss_i_tr = zeros(length(cSIGMAd_grid),1);
rss_i_dr = zeros(length(cSIGMAd_grid),1);


e_elb_tr = zeros(length(cSIGMAd_grid),1);
e_elb_dr = zeros(length(cSIGMAd_grid),1);

sim = 1000000;
burn = 1000;

for i = 1:length(cSIGMAd_grid)
    cSIGMAd = cSIGMAd_grid(i);
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd];
    
    %% Define i.i.d Discretizing State-Space
    [R.del_grid,R.s]=iid(cSIGMAd,nstate);
    
    %% Target Regime
    % Get initial matrix and solution vector
    [A_tr, b_tr, out_tr] = eqmmat_tr(params,R,nstate);

    y_tr = out_tr(1:nstate);
    pi_tr = out_tr(nstate+1:2*nstate);
    i_tr = out_tr(2*nstate+1:3*nstate);

    % Refine Matrix to account for ZLB
    if sum(i_tr >= 0) == nstate
        converged = 1;
    else
        converged = 0;
    end

    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_tr(params,R,nstate,A_tr,b_tr);

        i_tr = out_up(2*nstate+1:3*nstate);

        if sum(i_tr >= 0) == nstate
            converged =1;

            y_tr = out_up(1:nstate);
            pi_tr = out_up(nstate+1:2*nstate);
            i_tr = out_up(2*nstate+1:3*nstate);
            for k = 1:size(i_tr,1)
                if abs(i_tr(k)) < 1e-10
                    i_tr(k) = 0;
                end
            end
        end

        A_tr = A_up;
        b_tr = b_up;
    end
    
    % Get RSS
    rss_y_tr(i) = y_tr(rss_inx);
    rss_pi_tr(i) = pi_tr(rss_inx);
    rss_i_tr(i) = i_tr(rss_inx);
    
    % Get Expected Values
    pf_tr = [y_tr, pi_tr, i_tr];
    [e_tr,e_elb_t] = exp_val(params,pf_tr,R,sim,burn,1);

    e_y_tr(i) = e_tr(1);%exp_val(params,y_tr,R,sim,burn);
    e_pi_tr(i) = e_tr(2);%exp_val(params,pi_tr,R,sim,burn);
    e_i_tr(i) = e_tr(3);%e_elb_tr(i)] = exp_val(params,i_tr,R,sim,burn,1);
    e_elb_tr(i) = e_elb_t;


    
    %% Deflationary Regime
    % Get initial matrix and solution vector
    [A_dr, b_dr, out_dr] = eqmmat_dr(params,R,nstate);

    y_dr = out_dr(1:nstate);
    pi_dr = out_dr(nstate+1:2*nstate);
    i_dr = out_dr(2*nstate+1:3*nstate);

    i_real_dr = cRstar + cPHIpi*pi_dr;

    % Refine Matrix to account for ZLB
    if sum(i_dr > 0) == sum(i_real_dr > 0)
        converged = 1;
    else
        converged = 0;
    end
    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_dr(params,R,nstate,A_dr,b_dr);

        pi_dr = out_up(nstate+1:2*nstate);
        i_dr = out_up(2*nstate+1:3*nstate);
        i_real_dr = cRstar + cPHIpi*pi_dr;


        if sum(i_dr > 0) == sum(i_real_dr > 0)
            converged =1;

            y_dr = out_up(1:nstate);
            pi_dr = out_up(nstate+1:2*nstate);
            i_dr = out_up(2*nstate+1:3*nstate);
            for k = 1:size(i_dr,1)
                if abs(i_dr(k)) < 1e-10
                    i_dr(k) = 0;
                end
            end
        end

        A_dr = A_up;
        b_dr = b_up;
    end
    
    % Get RSS
    rss_y_dr(i) = y_dr(rss_inx);
    rss_pi_dr(i) = pi_dr(rss_inx);
    rss_i_dr(i) = i_dr(rss_inx);
    
    % Get Expected Values
    pf_dr = [y_dr, pi_dr, i_dr];
    [e_dr,e_elb_d] = exp_val(params,pf_dr,R,sim,burn,1);

    e_y_dr(i) = e_dr(1);%exp_val(params,y_dr,R,sim,burn);
    e_pi_dr(i) = e_dr(2);%exp_val(params,pi_dr,R,sim,burn);
    e_i_dr(i) = e_dr(3);%, e_elb_dr(i)] = exp_val(params,i_dr,R,sim,burn,1);
    e_elb_dr(i) = e_elb_d;
    
    disp(strcat('cSIGMAd = ',num2str(i),'of',num2str(length(cSIGMAd_grid))))
end

%% Plotting
X_tr = [100*rss_y_tr,400*rss_pi_tr,400*rss_i_tr,100*e_y_tr,400*e_pi_tr,400*e_i_tr,100*e_elb_tr];
X_dr = [100*rss_y_dr,400*rss_pi_dr,400*rss_i_dr,100*e_y_dr,400*e_pi_dr,400*e_i_dr,100*e_elb_dr];
ylims = [0 0.05; -1 0; 0 1;-0.01 0.02; -1 0.001; 0 1.001;0 100];

header = {'y_{RSS}','\pi_{RSS}','i_{RSS}','E[y]','E[\pi]','E[i]','Prob[ELB]'};

fig(1) = figure(1);
for i = 1:length(header);
    subplot(3,3,i)
    box on
    grid on
    hold on
    if i == length(header);
        h(1) = plot(cSIGMAd_grid,X_tr(:,i),'Color','k','LineStyle','-','LineWidth',2);
        h(2) = plot(cSIGMAd_grid,X_dr(:,i),'Color','k','LineStyle','--','LineWidth',2);
    else
        plot(cSIGMAd_grid,X_tr(:,i),'Color','k','LineStyle','-','LineWidth',2)
        plot(cSIGMAd_grid,X_dr(:,i),'Color','k','LineStyle','--','LineWidth',2)
    end
    
    xlabel('\sigma_{\epsilon}','FontSize',15)
    title(header{i},'FontSize',15,'FontWeight','Normal')
    set(gca,'Xlim',[cSIGMAd_grid(2) cSIGMAd_grid(end)],'Ylim',ylims(i,:),'FontSize',15)
end
L = legend([h(1) h(2)], 'Target Equilibrium', 'Deflationary Equilibrium');
set(L,'Position',[0.475 0.2 0.1 0.1],'FontSize',15)

savedir = cd;
if ispc 
    savedir = fullfile(savedir, '..\..\..');
    savedir = strcat(savedir,'\Final\');
else
    savedir = fullfile(savedir, '../../..');
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','stylized_iid_moments_cPHIpi_l.eps');
print(fig(1),'-depsc',strcat(savedir,'stylized_iid_moments_cPHIpi_l.eps'));


