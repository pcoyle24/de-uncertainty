% --------------------------------------------------------------------------
% File Name: script.m
% Author: Philip Coyle
% Date Created: 02/06/2018
% cd P:\RA_Work\Taisuke_Nakata\Zero_Lower_Bound\DeflationaryRegime\Cheb_play\sunspot\1d\ZLB
% --------------------------------------------------------------------------

% Projection Method (Chebyshev basis):
% Canonical New Keynesian Model (Rotemberg Pricing) with AR1 Preference
% Shock


clear all;
close all;
clc;
tstart = tic;                           % Job timer start

% dbstop in eqm at 37 if i>=6
% dbstop if error

%--------------------------------------------------------------------------
% Initialize Policy Functions
%--------------------------------------------------------------------------
% Load parameters, steady state and grids
P = parameters;

% Specify grid options
bound = P.sigma/((1-P.rho)^0.5);
O.delbound = [1-4*bound 1+4*bound];
O.del_pts = 11;
O.sunbound = [P.Pd, P.Ps];
O.sun_pts = 2;
O.e_pts = 10;

P.range = 0*(chebspace(O.delbound(2),O.delbound(1),O.del_pts)'-1);

% Chebyshev polynomial order in each dimension (must <= pts in grid)
O.n1 = 4;

% Conditional Probability for Sun Spot Shock
TransMat = [P.Ps, 1-P.Ps;
            1-P.Pd, P.Pd]; 

% Load discrsetized state space
G = grids_cheb(O,P);

% Initialize Chebyshev polynomial parameters
C = chebpoly(G,O);

global pi_yesterday
pi_yesterday       = 1;
for i = 1:length(P.pi_targ)
    S = steadystate(P,i);
    
    % Retrieve initial policy functions
    run_first_time = 1; %1 = True; 0 = False
    pf = guess_nzlb(P,S,G,O,C,run_first_time); 

    %----------------------------------------------------------------------
    % Initialize algorithm
    %----------------------------------------------------------------------
    % Calculate initial Chebyshev coefficients (Standard Regime, Non-ZLB)          
    C.Ac_s = chebweights11(O.n1,pf.c_s,C.T,C.P,C.X);
    C.Ainf_s = chebweights11(O.n1,pf.inf_s,C.T,C.P,C.X);

    % Calculate initial Chebyshev coefficients (Stanard Regime, ZLB)          
    C.Ac_zlb_s = chebweights11(O.n1,pf.c_zlb_s,C.T,C.P,C.X);
    C.Ainf_zlb_s = chebweights11(O.n1,pf.inf_zlb_s,C.T,C.P,C.X);

    % Calculate initial Chebyshev coefficients (Deflationary Regime, Non-ZLB)          
    C.Ac_d = chebweights11(O.n1,pf.c_d,C.T,C.P,C.X);
    C.Ainf_d = chebweights11(O.n1,pf.inf_d,C.T,C.P,C.X);

    % Calculate initial Chebyshev coefficients (Deflationary Regime, ZLB)          
    C.Ac_zlb_d = chebweights11(O.n1,pf.c_zlb_d,C.T,C.P,C.X);
    C.Ainf_zlb_d = chebweights11(O.n1,pf.inf_zlb_d,C.T,C.P,C.X);


    theta_old = [C.Ac_s,C.Ainf_s,C.Ac_zlb_s,C.Ainf_zlb_s,C.Ac_d,C.Ainf_d,C.Ac_zlb_d,C.Ainf_zlb_d];
    optionsj = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',100,'MaxIter',100,'TolFun',1e-10,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',true);

    
    func    = @(weights) eqm_jacob(weights,P,S,G,O,C,TransMat);
    display(char('Calculating the Policy Functions for Economy with ZLB'))
    tic
    [theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func,theta_old,[],[],optionsj);
    toc
    
%     options = optimoptions('lsqnonlin','display','iter','MaxFunEvals',100000000,'MaxIter',100000000,'TolFun',1e-9,'TolX',1e-9,'FunValCheck','On');
% 
%     func    = @(weights) eqm(weights,P,S,G,O,C,TransMat);
%     display(char('Calculating the Policy Functions for Economy with ZLB'))
%     [theta_new,out.norm,out.residual,out.exitflag,out.output,out.lambda,out.jacob] = lsqnonlin(func,theta_old,[],[],options);

   
    % Weights (Standard Regime)
    C.Ac_s        = theta_new(1:O.n1+1,1);
    C.Ainf_s      = theta_new(1:O.n1+1,2);
    C.Ac_zlb_s    = theta_new(1:O.n1+1,3);
    C.Ainf_zlb_s  = theta_new(1:O.n1+1,4);
    % Weights (Deflationary Regime)
    C.Ac_d        = theta_new(1:O.n1+1,5);
    C.Ainf_d      = theta_new(1:O.n1+1,6);
    C.Ac_zlb_d    = theta_new(1:O.n1+1,7);
    C.Ainf_zlb_d  = theta_new(1:O.n1+1,8);

    % Build out Policy functions
    % Standard Regime, NZLB
    pf.c_s = allcheb111(O.delbound,G.del,C.Ac_s,C.T,C.P);    
    pf.inf_s = allcheb111(O.delbound,G.del,C.Ainf_s,C.T,C.P); 
    pf.pitilde_s  = pf.inf_s/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n_s        = (pf.c_s./(1-(P.varphi/2).*(pf.pitilde_s-1).^2));
    pf.y_s        = pf.n_s;
    pf.w_s        = pf.n_s.^P.chin.*pf.c_s.^P.chic;
    pf.r_s        = S.pi_targ/P.beta*((pf.inf_s./S.pi_targ).^(P.phi_pi).*(pf.y_s./S.y).^(P.phi_y));
    
    % Standard Regime, ZLB
    pf.c_zlb_s = allcheb111(O.delbound,G.del,C.Ac_zlb_s,C.T,C.P);    
    pf.inf_zlb_s = allcheb111(O.delbound,G.del,C.Ainf_zlb_s,C.T,C.P); 
    pf.pitilde_zlb_s  = pf.inf_zlb_s/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n_zlb_s        = (pf.c_zlb_s./(1-(P.varphi/2).*(pf.pitilde_zlb_s-1).^2));
    pf.y_zlb_s        = pf.n_zlb_s;
    pf.w_zlb_s        = pf.n_zlb_s.^P.chin.*pf.c_zlb_s.^P.chic;
    pf.r_zlb_s = ones(length(pf.r_s),1);
    
    % Deflationary Regime, NZLB
    pf.c_d = allcheb111(O.delbound,G.del,C.Ac_d,C.T,C.P);    
    pf.inf_d = allcheb111(O.delbound,G.del,C.Ainf_d,C.T,C.P); 
    pf.pitilde_d  = pf.inf_d/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n_d        = (pf.c_d./(1-(P.varphi/2).*(pf.pitilde_d-1).^2));
    pf.y_d        = pf.n_d;
    pf.w_d        = pf.n_d.^P.chin.*pf.c_d.^P.chic;
    pf.r_d        = S.pi_targ/P.beta*((pf.inf_d./S.pi_targ).^(P.phi_pi).*(pf.y_d./S.y).^(P.phi_y)); 
    
    % Deflationary Regime, ZLB
    pf.c_zlb_d = allcheb111(O.delbound,G.del,C.Ac_zlb_d,C.T,C.P);    
    pf.inf_zlb_d = allcheb111(O.delbound,G.del,C.Ainf_zlb_d,C.T,C.P); 
    pf.pitilde_zlb_d  = pf.inf_zlb_d/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n_zlb_d        = (pf.c_zlb_d./(1-(P.varphi/2).*(pf.pitilde_zlb_d-1).^2));
    pf.y_zlb_d        = pf.n_zlb_d;
    pf.w_zlb_d        = pf.n_zlb_d.^P.chin.*pf.c_zlb_d.^P.chic;
    pf.r_zlb_d = ones(length(pf.r_d),1);

    display(char('Policy Functions Calculated'))
    
    save('pf_sun.mat','pf','O','P','S','G','C');
end



%--------------------------------------------------------------------------
% Construct CF Policy Functions
%--------------------------------------------------------------------------

pf.cf.c_s = zeros(G.griddim);
pf.cf.inf_s = zeros(G.griddim);
pf.cf.pitilde_s = zeros(G.griddim);
pf.cf.n_s = zeros(G.griddim);
pf.cf.y_s = zeros(G.griddim);
pf.cf.w_s = zeros(G.griddim);
pf.cf.r_s = zeros(G.griddim);

pf.cf.c_d = zeros(G.griddim);
pf.cf.inf_d = zeros(G.griddim);
pf.cf.pitilde_d = zeros(G.griddim);
pf.cf.n_d = zeros(G.griddim);
pf.cf.y_d = zeros(G.griddim);
pf.cf.w_d = zeros(G.griddim);
pf.cf.r_d = zeros(G.griddim);

for i = 1:G.nodes
    if pf.r_s(i) >= 1 
        pf.cf.c_s(i) = pf.c_s(i);
        pf.cf.inf_s(i) = pf.inf_s(i);
        pf.cf.pitilde_s(i) = pf.pitilde_s(i);
        pf.cf.n_s(i) = pf.n_s(i);
        pf.cf.y_s(i) = pf.y_s(i);
        pf.cf.w_s(i) = pf.w_s(i);
        pf.cf.r_s(i) = pf.r_s(i);
    else
        pf.cf.c_s(i) = pf.c_zlb_s(i);
        pf.cf.inf_s(i) = pf.inf_zlb_s(i);
        pf.cf.pitilde_s(i) = pf.pitilde_zlb_s(i);
        pf.cf.n_s(i) = pf.n_zlb_s(i);
        pf.cf.y_s(i) = pf.y_zlb_s(i);
        pf.cf.w_s(i) = pf.w_zlb_s(i);
        pf.cf.r_s(i) = pf.r_zlb_s(i);
    end
    if pf.r_d(i) >= 1 
        pf.cf.c_d(i) = pf.c_d(i);
        pf.cf.inf_d(i) = pf.inf_d(i);
        pf.cf.pitilde_d(i) = pf.pitilde_d(i);
        pf.cf.n_d(i) = pf.n_d(i);
        pf.cf.y_d(i) = pf.y_d(i);
        pf.cf.w_d(i) = pf.w_d(i);
        pf.cf.r_d(i) = pf.r_d(i);
    else
        pf.cf.c_d(i) = pf.c_zlb_d(i);
        pf.cf.inf_d(i) = pf.inf_zlb_d(i);
        pf.cf.pitilde_d(i) = pf.pitilde_zlb_d(i);
        pf.cf.n_d(i) = pf.n_zlb_d(i);
        pf.cf.y_d(i) = pf.y_zlb_d(i);
        pf.cf.w_d(i) = pf.w_zlb_d(i);
        pf.cf.r_d(i) = pf.r_zlb_d(i);
    end
end 


%--------------------------------------------------------------------------
% Plot Policy Functions
%-------------------------------------------------------------------------- 
zlb_s = find(pf.r_s < 1,1); %cell where ZLB binds for increasing \delta levels

Xs = [(pf.cf.c_s-S.c_s)./S.c_s*100, 400*(pf.cf.inf_s-1), (pf.cf.y_s-S.y_s)./S.y_s*100, 400*(pf.cf.r_s-1)];
disp(400*(pf.cf.r_s-1));
Xd = [(pf.cf.c_d-S.c_s)./S.c_s*100, 400*(pf.cf.inf_d-1), (pf.cf.y_d-S.y_s)./S.y_s*100, 400*(pf.cf.r_d-1)];
title_fig_s = {'Consumption (S)', 'Inflation Rate (S)', 'Output (S)', 'Interest Rate (S)'};
title_fig_d = {'Consumption (D)', 'Inflation Rate (D)', 'Output (D)', 'Interest Rate (D)'};

ss = [(S.c_d-S.c_s)/S.c_s*100 (S.c_d-S.c_s)/S.c_s*100; 0 0; (S.y_d-S.y_s)/S.y_s*100 (S.y_d-S.y_s)/S.y_s*100; 0 0];

% ylim_s = [-.075 .05;-10 12;-.05 .15;0 20];
% ylim_d = [-.075 .02;-10 -0;-.04 .01;-1 1];

fig(1) = figure(1);
for i = 1:length(title_fig_s)
    subplot(2,2,i)
    box on
    hold on
    grid on
    plot(G.del,Xs(:,i),'k','LineWidth',2)
    xlabel('\delta','FontSize',16)
    title(title_fig_s{i},'FontSize',16)
%    set(gca,'XLim',[1-4*bound 1+4*bound],'YLim',ylim_s(i,:),'FontSize',16)
    set(gca,'XLim',[1-4*P.bound 1+4*P.bound],'FontSize',16)
%    line([G.del(zlb_s,1) G.del(zlb_s,1)],get(gca,'YLim'),'Color','b','LineStyle','--','LineWidth',1.5)
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','Pfs_SunSpotEconomy_S.eps');

fig(2) = figure(2);
for i = 1:length(title_fig_d)
    subplot(2,2,i)
    box on
    hold on
    grid on
    plot(G.del,Xd(:,i),'k','LineWidth',2)
    xlabel('\delta','FontSize',16)
    title(title_fig_d{i},'FontSize',16)
%    set(gca,'XLim',O.delbound,'YLim',ylim_d(i,:),'FontSize',16)
    set(gca,'XLim',O.delbound,'FontSize',16)
%    line(get(gca,'XLim'),ss(i,:),'Color','b','LineStyle','--','LineWidth',1.5)
end

set(fig(2),'PaperOrientation','Landscape');
set(fig(2),'PaperPosition',[0 0 11 8.5]);
% print(fig(2),'-depsc','Pfs_SunSpotEconomy_D.eps');