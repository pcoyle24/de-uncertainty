% --------------------------------------------------------------------------
% File Name: stdev_pfs.m
% Author: Philip Coyle
% Date Created: 12/11/2018
% 
% Add Appropriate Paths 
% codetoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/Code_Toolbox'; 
% mextoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/mex_functions';
% addpath(genpath(codetoolbox));
% addpath(genpath(mextoolbox)); 
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/AR1/1d/sunspot/Datagen
% stdev_pfs
% --------------------------------------------------------------------------

clear all
close all 
clc

%% Housekeeping
if ~ispc
    % start parallel pool
    gcp;
end

% Set Seed
rng(15)

% Add appropriate Paths
addpath savedata/

%% Calcuate moments 
% Load parameters
P = parameters;
if ~ispc
    % Allocate space
    E_c = zeros(length(P.sigma_grid),2);
    E_pi = zeros(length(P.sigma_grid),2);
    E_r = zeros(length(P.sigma_grid),2);

    elb_prob = zeros(length(P.sigma_grid),2);

    sig_c = zeros(length(P.sigma_grid),2);
    sig_pi = zeros(length(P.sigma_grid),2);
    sig_r = zeros(length(P.sigma_grid),2);

    rss_c_s = zeros(length(P.sigma_grid),2);
    rss_inf_s = zeros(length(P.sigma_grid),2);
    rss_r_s = zeros(length(P.sigma_grid),2);
    rss_c_d = zeros(length(P.sigma_grid),2);
    rss_inf_d = zeros(length(P.sigma_grid),2);
    rss_r_d = zeros(length(P.sigma_grid),2);

    rss_gap_c_s = zeros(length(P.sigma_grid),2);
    rss_gap_inf_s = zeros(length(P.sigma_grid),2);
    rss_gap_r_s = zeros(length(P.sigma_grid),2);
    rss_gap_c_d = zeros(length(P.sigma_grid),2);
    rss_gap_inf_d = zeros(length(P.sigma_grid),2);
    rss_gap_r_d = zeros(length(P.sigma_grid),2);

    for k = 1:2
        P.Ps = P.Ps_grid(k);
        P.Pd = P.Pd_grid(k);
        if k == 1
            disp('Let p_D be an absorbing state')
        else
            disp('Let p_T be an absorbing state')
        end

        for i = 1:length(P.sigma_grid)
            P.sigma = P.sigma_grid(i);
            P.bound = (P.sigma^2/(1-P.rho^2))^0.5;
            shock = 4*P.bound;

            load(strcat('ChebPFs_cSIGMA',num2str(100*P.sigma),'_ps',num2str(P.Ps),'_pd',num2str(P.Pd),'.mat'));

            [E, sig, Pelb] = sim_pfs_zlb(P,O,C,S);
            % Expected Value
            E_c(i,k) = E.c;
            E_pi(i,k) = E.pi;
            E_r(i,k) = E.r;

            % ELB Probability
            elb_prob(i,k) = 100*Pelb;

            % Standard Deviation
            sig_c(i,k) = sig.c;
            sig_pi(i,k) = sig.pi;
            sig_r(i,k) = sig.r;

            disp(strcat('Unconditional Expected Values and Standard Deviations calculated for cSIGMA =',num2str(P.sigma)))
            disp(' ')

            % Record RSS and RSS Gap values
            rss_c_s(i,k) = rss.c_s;
            rss_inf_s(i,k) = rss.inf_s;
            rss_r_s(i,k) = rss.r_s;

            rss_c_d(i,k) = rss.c_d;
            rss_inf_d(i,k) = rss.inf_d;
            rss_r_d(i,k) = rss.r_d;

            rss_gap_c_s(i,k) = rss.gap.c_s;
            rss_gap_inf_s(i,k) = rss.gap.inf_s;
            rss_gap_r_s(i,k) = rss.gap.r_s;

            rss_gap_c_d(i,k) = rss.gap.c_d;
            rss_gap_inf_d(i,k) = rss.gap.inf_d;
            rss_gap_r_d(i,k) = rss.gap.r_d;

        end
    end

    save('momentsdata.mat','rss_c_s','rss_c_d','E_c','sig_c','rss_inf_s','rss_inf_d','E_pi','sig_pi','rss_r_s','rss_r_d','E_r','sig_r','elb_prob');
else
    load momentsdata.mat

    tr_mat = round([rss_c_s(:,2),E_c(:,2),sig_c(:,2),rss_inf_s(:,2),E_pi(:,2),sig_pi(:,2),rss_r_s(:,2),E_r(:,2),sig_r(:,2),elb_prob(:,2)],2);
    tr_mat = [P.sigma_grid', tr_mat];

    dr_mat = round([rss_c_d(:,1),E_c(:,1),sig_c(:,1),rss_inf_d(:,1),E_pi(:,1),sig_pi(:,1),rss_r_d(:,1),E_r(:,1),sig_r(:,1),elb_prob(:,1)],2);
    dr_mat = [P.sigma_grid', dr_mat];

    title = {'cSIGMA','C_{RSS}','E{C}','\sigma_{C}','\Pi_{RSS}','E{\Pi}','\sigma_{\Pi}','R_{RSS}','E{R}','\sigma_{R}','P{ELB}'};
    regime = {'';'Target Regime';'';'';'* Largest shock such that RSS_{R}^{D} = 0';'';'';'Deflationary Regime';'';'';'* Largest shock such that RSS_{R}^{D} = 0';'';'';};
    blankspace = {'','','','','','','','','','',''};

    out_mat = [title;...
          blankspace;...
          num2cell(tr_mat);...
          blankspace;...
    	  num2cell(dr_mat)];
    out_mat = [regime,out_mat];
      
    xlswrite('momentsdata.xlsx',out_mat);

end