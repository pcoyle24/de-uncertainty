% --------------------------------------------------------------------------
% File Name: plot_ti_pfs.m
% Author: Philip Coyle
% Date Created: 11/05/2018
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/AR1/1d/ZLB/Fortran

% plot_ti_pfs
% --------------------------------------------------------------------------

clear all
close all
clc

%% House Keeping
addpath('..')
P = parameters;
S = steadystate(P,1);

% bound = P.sigma/((1-P.rho)^0.5);
bound = (P.sigma^2/(1-P.rho^2))^0.5;
%% Load in Data
del_col = 1;
c_col = 2; 
pi_col = 3; 
n_col = 4; 
y_col = 5; 
w_col = 6; 
v_col = 7; 
r_col = 8; 

pfs_nzlb = dlmread('PFs_sp_d_nzlb.dat');
pfs_zlb = dlmread('PFs_sp_d_zlb.dat');

del_grid = pfs_nzlb(:,del_col);
pf_c  = pfs_nzlb(:,c_col);
pf_inf  = pfs_nzlb(:,pi_col);
pf_y  = pfs_nzlb(:,y_col);
pf_v  = pfs_nzlb(:,v_col);
pf_r  = pfs_nzlb(:,r_col);

pf_c_zlb  = pfs_zlb(:,c_col);
pf_inf_zlb  = pfs_zlb(:,pi_col);
pf_y_zlb  = pfs_zlb(:,y_col);
pf_v_zlb  = pfs_zlb(:,v_col);
pf_r_zlb  = pfs_zlb(:,r_col);

%% CF Construction
del_grid_fine = linspace(del_grid(1), del_grid(end), 1001)';
pf_cf_c = zeros(length(del_grid_fine),1);
pf_cf_inf = zeros(length(del_grid_fine),1);
pf_cf_v = zeros(length(del_grid_fine),1);
pf_cf_n = zeros(length(del_grid_fine),1);
pf_cf_y = zeros(length(del_grid_fine),1);
pf_cf_w = zeros(length(del_grid_fine),1);
pf_cf_v = zeros(length(del_grid_fine),1);
pf_cf_r = zeros(length(del_grid_fine),1);


for i = 1:length(del_grid_fine)
    zlb_check = allinterp1(del_grid, del_grid_fine(i),pf_r);
    if zlb_check >= 0 
        pf_cf_c(i) = allinterp1(del_grid, del_grid_fine(i),pf_c);
        pf_cf_inf(i) = allinterp1(del_grid, del_grid_fine(i),pf_inf);
        pf_cf_y(i) = allinterp1(del_grid, del_grid_fine(i),pf_y);
        pf_cf_v(i) = allinterp1(del_grid, del_grid_fine(i),pf_v);
        pf_cf_r(i) = allinterp1(del_grid, del_grid_fine(i),pf_r);
    else
        pf_cf_c(i) = allinterp1(del_grid, del_grid_fine(i),pf_c_zlb);
        pf_cf_inf(i) = allinterp1(del_grid, del_grid_fine(i),pf_inf_zlb);
        pf_cf_y(i) = allinterp1(del_grid, del_grid_fine(i),pf_y_zlb);
        pf_cf_v(i) = allinterp1(del_grid, del_grid_fine(i),pf_v_zlb);
        pf_cf_r(i) = allinterp1(del_grid, del_grid_fine(i),pf_r_zlb);
    end
end

%% Plot Policy Functions

zlb = find(pf_cf_r == 0,1); %cell where ZLB binds for increasing \delta levels

% X = [100*(pf_cf_c-S.c)./S.c, pf_cf_inf, 100*(pf_cf_y-S.y)./S.y, pf_cf_r];
X = [100*(pf_cf_c-S.c)./S.c, pf_cf_inf, pf_cf_r];
% disp(400*(pf_r-1));
title_fig = {'Consumption', 'Inflation Rate', 'Interest Rate'};

% ylim = [-15 10;-8 10;-250 -190;0 15];

fig(1) = figure(1);
for i = 1:length(title_fig)
    subplot(2,3,i)
    box on
    hold on
    grid on
    plot(del_grid_fine,X(:,i),'k','LineWidth',2)
    xlabel('\delta','FontSize',16)
    title(title_fig{i},'FontSize',16)
%     set(gca,'XLim',[1-4*bound 1+4*bound],'Ylim',ylim(i,:),'FontSize',16)
    set(gca,'XLim',[1-4*bound 1+4*bound],'FontSize',16)
    line([del_grid_fine(zlb,1) del_grid_fine(zlb,1)],get(gca,'YLim'),'Color','b','LineStyle','--','LineWidth',1.5)
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','Pfs_TargetRegime_RSSZLBBind.eps');
