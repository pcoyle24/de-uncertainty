%--------------------------------------------------------------------------
% File Name: stylized_rouwenhorst_moments_cPHIpi_;.m
% Author: Philip Coyle
% Date Created: 01/18/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/ar1/nstate/moments
% matlab -nodesktop -nosplash -r stylized_rouwenhorst_moments_cPHIpi_l
%--------------------------------------------------------------------------

clear all
close all
clc

load('cSIGMAd_max.mat')
%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 2;
cRstar = 1/400;
cRHO = 0.8; 
cPItarg= 2/400;
cSIGMAd_grid = linspace(0,cSIGMAd_max.l-eps,51);

numpts = 21;

%% Housekeeping
mu_y_tr = zeros(length(cSIGMAd_grid),1);
mu_y_dr = zeros(length(cSIGMAd_grid),1);
sig2_y_tr = zeros(length(cSIGMAd_grid),1);
sig2_y_dr = zeros(length(cSIGMAd_grid),1);

rss_y_tr = zeros(length(cSIGMAd_grid),1);
rss_y_dr = zeros(length(cSIGMAd_grid),1);

mu_pi_tr = zeros(length(cSIGMAd_grid),1);
mu_pi_dr = zeros(length(cSIGMAd_grid),1);
sig2_pi_tr = zeros(length(cSIGMAd_grid),1);
sig2_pi_dr = zeros(length(cSIGMAd_grid),1);

rss_pi_tr = zeros(length(cSIGMAd_grid),1);
rss_pi_dr = zeros(length(cSIGMAd_grid),1);

mu_i_tr = zeros(length(cSIGMAd_grid),1);
mu_i_dr = zeros(length(cSIGMAd_grid),1);
sig2_i_tr = zeros(length(cSIGMAd_grid),1);
sig2_i_dr = zeros(length(cSIGMAd_grid),1);

rss_i_tr = zeros(length(cSIGMAd_grid),1);
rss_i_dr = zeros(length(cSIGMAd_grid),1);


mu_elb_tr = zeros(length(cSIGMAd_grid),1);
mu_elb_dr = zeros(length(cSIGMAd_grid),1);
sig2_elb_tr = zeros(length(cSIGMAd_grid),1);
sig2_elb_dr = zeros(length(cSIGMAd_grid),1);

sim = 1000000;
burn = 1000;

tic;
for i = 1:length(cSIGMAd_grid)
    cSIGMAd = cSIGMAd_grid(i);
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd;cPItarg];


    %% Define Rouwenhorst Discretizing State-Space
    [R.del_grid,R.trans_mat,R.s]=rouwenhorst(cRHO,cSIGMAd,numpts);
    if i == 1;
        rss_inx = (numpts+1)/2;
    else
        rss_inx = find(abs(R.del_grid - 0) < 1e-10);
    end

    %% Target Regime
    % Get initial matrix and solution vector
    [A_tr, b_tr, out_tr] = eqmmat_tr(params,R,numpts);

    y_tr = out_tr(1:numpts);
    pi_tr = out_tr(numpts+1:2*numpts);
    i_tr = out_tr(2*numpts+1:3*numpts);

    % Refine Matrix to account for ZLB
    if sum(i_tr >= 0) == numpts
        converged = 1;
    else
        converged = 0;
    end

    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_tr(params,R,numpts,A_tr,b_tr);

        i_tr = out_up(2*numpts+1:3*numpts);

        if sum(i_tr >= 0) == numpts
            converged =1;

            y_tr = out_up(1:numpts);
            pi_tr = out_up(numpts+1:2*numpts);
            i_tr = out_up(2*numpts+1:3*numpts);

            % Check equilibirum existence
            i_real_tr = cRstar + cPItarg + cPHIpi*(pi_tr - cPItarg);
            inx = find(i_real_tr>0);
            if sum(abs(i_tr(inx) - i_real_tr(inx)) < eps) ~= length(inx)
                error('The equilibrium does not exist. Reduce shock size.')
            end
        end

        A_tr = A_up;
        b_tr = b_up;
    end
    
    % Get RSS
    rss_y_tr(i) = y_tr(rss_inx);
    rss_pi_tr(i) = pi_tr(rss_inx);
    rss_i_tr(i) = i_tr(rss_inx);
    
    % Simulate Data
    y_tr_sim = sim_data(params,y_tr,R,sim,burn);
    pi_tr_sim = sim_data(params,pi_tr,R,sim,burn);
    [i_tr_sim, elb_tr_sim] = sim_data(params,i_tr,R,sim,burn,1);
    
    % Get Expected Values
    mu_y_tr(i) = mean(y_tr_sim);
    mu_pi_tr(i) = mean(pi_tr_sim);
    mu_i_tr(i) = mean(i_tr_sim);
    mu_elb_tr(i) = mean(elb_tr_sim);
    
    % Get Variances
    sig2_y_tr(i) = var(y_tr_sim);
    sig2_pi_tr(i) = var(pi_tr_sim);
    sig2_i_tr(i) = var(i_tr_sim);
    sig2_elb_tr(i) = var(elb_tr_sim);

%     % Get Expected Values
%     mu_y_tr(i) = exp_val(params,y_tr,R,sim,burn);
%     mu_pi_tr(i) = exp_val(params,pi_tr,R,sim,burn);
%     [mu_i_tr(i), mu_elb_tr(i)] = exp_val(params,i_tr,R,sim,burn,1);


    
    %% Deflationary Regime
    % Get initial matrix and solution vector
    [A_dr, b_dr, out_dr] = eqmmat_dr(params,R,numpts);

    y_dr = out_dr(1:numpts);
    pi_dr = out_dr(numpts+1:2*numpts);
    i_dr = out_dr(2*numpts+1:3*numpts);

    i_real_dr = cRstar + cPItarg + cPHIpi*(pi_dr - cPItarg);

    % Refine Matrix to account for ZLB
    if sum(i_dr > 0) == sum(i_real_dr > 0)
        converged = 1;
    else
        converged = 0;
    end
    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_dr(params,R,numpts,A_dr,b_dr);

        pi_dr = out_up(numpts+1:2*numpts);
        i_dr = out_up(2*numpts+1:3*numpts);
        i_real_dr = cRstar + cPItarg + cPHIpi*(pi_dr - cPItarg);


        if sum(i_dr > 0) == sum(i_real_dr > 0)
            converged =1;

            y_dr = out_up(1:numpts);
            pi_dr = out_up(numpts+1:2*numpts);
            i_dr = out_up(2*numpts+1:3*numpts);

            % Check equilibirum existence
            if sum(i_dr<0) > 0
                error('The equilibrium does not exist. Reduce shock size.');
            end
        end

        A_dr = A_up;
        b_dr = b_up;
    end
    
    % Get RSS
    rss_y_dr(i) = y_dr(rss_inx);
    rss_pi_dr(i) = pi_dr(rss_inx);
    rss_i_dr(i) = i_dr(rss_inx);
    
    % Simulate Data
    y_dr_sim = sim_data(params,y_dr,R,sim,burn);
    pi_dr_sim = sim_data(params,pi_dr,R,sim,burn);
    [i_dr_sim, elb_dr_sim] = sim_data(params,i_dr,R,sim,burn,1);
    
    % Get Expected Values
    mu_y_dr(i) = mean(y_dr_sim);
    mu_pi_dr(i) = mean(pi_dr_sim);
    mu_i_dr(i) = mean(i_dr_sim);
    mu_elb_dr(i) = mean(elb_dr_sim);
    
    % Get Variances
    sig2_y_dr(i) = var(y_dr_sim);
    sig2_pi_dr(i) = var(pi_dr_sim);
    sig2_i_dr(i) = var(i_dr_sim);
    sig2_elb_dr(i) = var(elb_dr_sim);

%     % Get Expected Values
%     mu_y_dr(i) = exp_val(params,y_dr,R,sim,burn);
%     mu_pi_dr(i) = exp_val(params,pi_dr,R,sim,burn);
%     [mu_i_dr(i), mu_elb_dr(i)] = exp_val(params,i_dr,R,sim,burn,1);
    
    disp(strcat('cSIGMAd = '," ",num2str(i)," ",'of'," ",num2str(length(cSIGMAd_grid))))
end
toc;

%% Plotting
%% Plotting
X_tr = [100*rss_y_tr,400*rss_pi_tr,400*rss_i_tr,...
        100*mu_y_tr,400*mu_pi_tr,400*mu_i_tr,...
        100^2*sig2_y_tr,400^2*sig2_pi_tr,400^2*sig2_i_tr,...
        100*mu_elb_tr];

X_dr = [100*rss_y_dr,400*rss_pi_dr,400*rss_i_dr,...
        100*mu_y_dr,400*mu_pi_dr,400*mu_i_dr,...
        100^2*sig2_y_dr,400^2*sig2_pi_dr,400^2*sig2_i_dr,...
        100*mu_elb_dr];

ylims = [0 0.4; -1 0; 0 1;-0.01 0.005; -1.01 0.1; 0 1.2;0 100];

header = {'y_{RSS}','\pi_{RSS}','i_{RSS}',...
          'E[y]','E[\pi]','E[i]',...
          'Var[y]','Var[\pi]','Var[i]',...
          'Prob[ELB]'};

fig(1) = figure(1);
for i = 1:length(header);
    subplot(4,3,i)
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
%     set(gca,'Xlim',[cSIGMAd_grid(2) cSIGMAd_grid(end)],'Ylim',ylims(i,:),'FontSize',15)
    set(gca,'Xlim',[cSIGMAd_grid(2) cSIGMAd_grid(end)],'FontSize',15)
end
L = legend([h(1) h(2)], 'Target Equilibrium', 'Deflationary Equilibrium');
set(L,'Position',[0.6 0.17 0.1 0.1],'FontSize',15)

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
print(fig(1),'-dpdf',strcat(savedir,'stylized_rouwenhorst_moments_cPHIpi_l.pdf'));