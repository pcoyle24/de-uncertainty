% --------------------------------------------------------------------------
% File Name: script_val_mex.m
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
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/AR1/1d/ZLB
% script_val_mex
% --------------------------------------------------------------------------

% Projection Method (Chebyshev basis):
% Canonical New Keynesian Model (Rotemberg Pricing) with AR1 Preference
% Shock


clear all;
close all;
clc;

% Directory to save data
savedir = pwd;
savedir = strcat(pwd,'/savedata/');

% Set Seed
rng(22)

%--------------------------------------------------------------------------
% Initialize Policy Functions
%--------------------------------------------------------------------------
% Load parameters and grids
P = parameters;

% Specify grid options
O.delbound = [1-4*P.bound 1+4*P.bound];
O.del_pts = 51;
O.e_pts = 10;

% Chebyshev polynomial order in each dimension (must <= pts in grid)
O.n1 = 4;

% Load discretized state space
G = grids_cheb(O,P);

% Initialize Chebyshev polynomial parameters
C = chebpoly(G,O);

global pi_yesterday
pi_yesterday       = 1;


for i = 1:length(P.pi_targ)
    S = steadystate(P,i);
    % Retrieve initial policy functions
    pf = guess_nzlb_val_mex(P,S,G,O,C);

    %--------------------------------------------------------------------------
    % Initialize algorithm
    %--------------------------------------------------------------------------

    % Calculate initial Chebyshev coefficients (Non-ZLB)          
    C.Ac = Fchebweights11(O.n1,O.del_pts,pf.c,C.T,C.P,C.X);
    C.Ainf = Fchebweights11(O.n1,O.del_pts,pf.inf,C.T,C.P,C.X);

    % Calculate initial Chebyshev coefficients (ZLB)          
    C.Ac_zlb = Fchebweights11(O.n1,O.del_pts,pf.c_zlb,C.T,C.P,C.X);
    C.Ainf_zlb = Fchebweights11(O.n1,O.del_pts,pf.inf_zlb,C.T,C.P,C.X);

    theta_old = [C.Ac,C.Ainf,C.Ac_zlb,C.Ainf_zlb];
    options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',300,'MaxIter',500,'TolFun',1e-12,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',true);
    func    = @(weights) eqm_mex(weights,P,S,G,O,C);

    display(char('Calculating the Policy Functions for Economy with ZLB'))

    [theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func,theta_old,[],[],options);

    C.Ac        = theta_new(1:O.n1+1,1);
    C.Ainf      = theta_new(1:O.n1+1,2);
    C.Ac_zlb    = theta_new(1:O.n1+1,3);
    C.Ainf_zlb  = theta_new(1:O.n1+1,4);

    pf.c        = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ac,C.max,C.T,C.P);    
    pf.inf      = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ainf,C.max,C.T,C.P); 
    pf.pitilde  = pf.inf/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n        = (pf.c./(1-(P.varphi/2).*(pf.pitilde-1).^2));
    pf.y        = pf.n;
    pf.w        = pf.n.^P.chin.*pf.c.^P.chic;
    pf.r        = S.pi_targ/P.beta*((pf.inf./S.pi_targ).^(P.phi_pi).*(pf.y./S.y).^(P.phi_y));

    pf.c_zlb        = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ac_zlb,C.max,C.T,C.P);    
    pf.inf_zlb      = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ainf_zlb,C.max,C.T,C.P);
    pf.pitilde_zlb  = pf.inf_zlb/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n_zlb        = (pf.c_zlb./(1-(P.varphi/2).*(pf.pitilde_zlb-1).^2));
    pf.y_zlb        = pf.n_zlb;
    pf.w_zlb        = pf.n_zlb.^P.chin.*pf.c_zlb.^P.chic;  
    pf.r_zlb        = ones(length(pf.r),1);

    C.Ar = Fchebweights11(O.n1,O.del_pts,pf.r,C.T,C.P,C.X);

    display(char('Policy Functions Calculated'))

    display(char('Calculating Value Function'))

    %----------------------------------------------------------------------
    % Initialize algorithm (Value Function Iteration) 
    %----------------------------------------------------------------------
    % Calculate initial Chebyshev coefficients (Non-ZLB)          
    C.Av = Fchebweights11(O.n1,O.del_pts,pf.v,C.T,C.P,C.X);
    % Calculate initial Chebyshev coefficients (ZLB)          
    C.Av_zlb = Fchebweights11(O.n1,O.del_pts,pf.v_zlb,C.T,C.P,C.X);

    theta_old_val = [C.Av,C.Av_zlb];

    % Preallocate arrays to store policy function updates    
    func_val    = @(weights) eqm_val(weights,P,S,G,O,C,pf);    
    options_val = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',300,'MaxIter',500,'TolFun',1e-12,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',false);

    [theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func_val,theta_old_val,[],[],options_val);

    % Map Coefficients 
    C.Av     = theta_new(1:O.n1+1,1);
    C.Av_zlb = theta_new(1:O.n1+1,2);

    pf.v =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av,C.max,C.T,C.P);    
    pf.v_zlb =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av_zlb,C.max,C.T,C.P); 
    
%     % Save coefficients
%     save(strcat(savedir,'ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'.mat'),'C','O','pf','S');

    %% ------------------------------------------------------------------------
    %  Construct CF Policy Functions
    %  ------------------------------------------------------------------------
    O.del_fine = 1001;
    G.delfine = chebspace(O.delbound(1),O.delbound(2),O.del_fine)';
    G.gridfine = size(ndgrid(G.delfine));


    pf_c_fine        = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ac,C.max,C.T,C.P);
    pf_inf_fine      = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ainf,C.max,C.T,C.P);
    pf_pitilde_fine  = pf_inf_fine/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf_n_fine        = (pf_c_fine./(1-(P.varphi/2).*(pf_pitilde_fine-1).^2));
    pf_y_fine        = pf_n_fine;
    pf_w_fine        = pf_n_fine.^P.chin.*pf_c_fine.^P.chic;
    pf_r_fine        = S.pi_targ/P.beta*((pf_inf_fine./S.pi_targ).^(P.phi_pi).*(pf_y_fine./S.y).^(P.phi_y));
    pf_v_fine        = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Av,C.max,C.T,C.P);

    pf_c_zlb_fine        = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ac_zlb,C.max,C.T,C.P);
    pf_inf_zlb_fine      = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ainf_zlb,C.max,C.T,C.P);
    pf_pitilde_zlb_fine  = pf_inf_zlb_fine/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf_n_zlb_fine        = (pf_c_zlb_fine./(1-(P.varphi/2).*(pf_pitilde_zlb_fine-1).^2));
    pf_y_zlb_fine        = pf_n_zlb_fine;
    pf_w_zlb_fine        = pf_n_zlb_fine.^P.chin.*pf_c_zlb_fine.^P.chic;
    pf_r_zlb_fine        = ones(length(pf_r_fine),1);
    pf_v_zlb_fine        = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Av_zlb,C.max,C.T,C.P);

    zlb = find(pf_r_fine < 1,1); %cell where ZLB binds for increasing \delta levels

    for j = 1:O.del_fine
        if pf_r_fine(j) >= 1 
            pf.cf.c(j) = pf_c_fine(j);
            pf.cf.inf(j) = pf_inf_fine(j);
            pf.cf.pitilde(j) = pf_pitilde_fine(j);
            pf.cf.n(j) = pf_n_fine(j);
            pf.cf.y(j) = pf_y_fine(j);
            pf.cf.w(j) = pf_w_fine(j);
            pf.cf.r(j) = pf_r_fine(j);
            pf.cf.v(j) = pf_v_fine(j);
        else
            pf.cf.c(j) = pf_c_zlb_fine(j);
            pf.cf.inf(j) = pf_inf_zlb_fine(j);
            pf.cf.pitilde(j) = pf_pitilde_zlb_fine(j);
            pf.cf.n(j) = pf_n_zlb_fine(j);
            pf.cf.y(j) = pf_y_zlb_fine(j);
            pf.cf.w(j) = pf_w_zlb_fine(j);
            pf.cf.r(j) = pf_r_zlb_fine(j);
            pf.cf.v(j) = pf_v_zlb_fine(j);
        end
    end

    %% ------------------------------------------------------------------------
    %  Plot Policy Functions
    %  ------------------------------------------------------------------------
    X = [100*(pf.cf.c'-S.c)./S.c, 400*(pf.cf.inf'-1), 400*(pf.cf.r'-1)];
    colors = {'k', 'b'};

    title_fig = {'Consumption', 'Inflation Rate', 'Interest Rate'};

    ylim = [-10 5;-4 6;0 12];

    fig(1) = figure(1);
    for k = 1:length(title_fig)
        subplot(2,3,k)
        box on
        hold on
        grid on
        if k == 1
            h(i) = plot(G.delfine,X(:,k),'Color',colors{i},'LineStyle','-','LineWidth',2);
        else
            plot(G.delfine,X(:,k),'Color',colors{i},'LineStyle','-','LineWidth',2);
        end
        xlabel('\delta','FontSize',16)
        title(title_fig{k},'FontSize',16,'FontWeight','Normal')
%         set(gca,'XLim',[1-4*P.bound 1+4*P.bound],'FontSize',16)
        set(gca,'XLim',[1-4*P.bound 1+4*P.bound],'YLim',ylim(k,:),'FontSize',16)
        line([G.delfine(zlb,1) G.delfine(zlb,1)],get(gca,'YLim'),'Color',colors{i},'LineStyle','--','LineWidth',1.5)
    end
end
% L = legend([h(1) h(2)],'\Pi^{targ} = 1.5%', '\Pi^{targ} = 2%');
% set(L,'Location','SouthWest','FontSize',15)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','Pfs.eps');
