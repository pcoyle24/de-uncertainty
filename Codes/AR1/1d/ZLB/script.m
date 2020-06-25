% --------------------------------------------------------------------------
% File Name: script.m
% Author: Philip Coyle
% Date Created: 02/06/2018
% Last Updated: 07/19/2018
% 
% cd P:\RA_Work\Taisuke_Nakata\Zero_Lower_Bound\DeflationaryRegime\Codes\StylizedModel\1d\ZLB
% script
% --------------------------------------------------------------------------

% Projection Method (Chebyshev basis):
% Canonical New Keynesian Model (Rotemberg Pricing) with AR1 Preference
% Shock


clear all;
close all;
clc;
tic

%--------------------------------------------------------------------------
% Initialize Policy Functions
%--------------------------------------------------------------------------
% Load parameters and grids
P = parameters1;

% Specify grid options
bound = P.sigma/((1-P.rho)^0.5);
O.delbound = [1-4*bound 1+4*bound];
O.del_pts = 11;
O.e_pts = 5;

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
    run_first_time = 1; %1 = True; 0 = False
    pf = guess_nzlb(P,S,G,O,C,run_first_time); 

    %--------------------------------------------------------------------------
    % Initialize algorithm
    %--------------------------------------------------------------------------
    % Calculate initial Chebyshev coefficients (Non-ZLB)          
    C.Ac = chebweights11(O.n1,pf.c,C.T,C.P,C.X);
    C.Ainf = chebweights11(O.n1,pf.inf,C.T,C.P,C.X);
    
    % Calculate initial Chebyshev coefficients (ZLB)          
    C.Ac_zlb = chebweights11(O.n1,pf.c_zlb,C.T,C.P,C.X);
    C.Ainf_zlb = chebweights11(O.n1,pf.inf_zlb,C.T,C.P,C.X);

    theta_old = [C.Ac,C.Ainf,C.Ac_zlb,C.Ainf_zlb];
    options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',true);
    func    = @(weights) eqm_jacob(weights,P,S,G,O,C);

    display(char('Calculating the Policy Functions for Economy with ZLB'))
    tic
    [theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func,theta_old,[],[],options);
    toc
    C.Ac        = theta_new(1:O.n1+1,1);
    C.Ainf      = theta_new(1:O.n1+1,2);
    C.Ac_zlb    = theta_new(1:O.n1+1,3);
    C.Ainf_zlb  = theta_new(1:O.n1+1,4);

    pf.c = allcheb111(O.delbound,G.del_grid,C.Ac,C.T,C.P);    
    pf.inf = allcheb111(O.delbound,G.del_grid,C.Ainf,C.T,C.P); 
    pf.pitilde  = pf.inf/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n        = (pf.c./(1-(P.varphi/2).*(pf.pitilde-1).^2));
    pf.y        = pf.n;
    pf.w        = pf.n.^P.chin.*pf.c.^P.chic;
    pf.r        = S.pi_targ/P.beta*((pf.inf./S.pi_targ).^(P.phi_pi).*(pf.y./S.y).^(P.phi_y));

    pf.c_zlb = allcheb111(O.delbound,G.del_grid,C.Ac_zlb,C.T,C.P);    
    pf.inf_zlb = allcheb111(O.delbound,G.del_grid,C.Ainf_zlb,C.T,C.P);
    pf.pitilde_zlb  = pf.inf_zlb/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
    pf.n_zlb        = (pf.c_zlb./(1-(P.varphi/2).*(pf.pitilde_zlb-1).^2));
    pf.y_zlb        = pf.n_zlb;
    pf.w_zlb        = pf.n_zlb.^P.chin.*pf.c_zlb.^P.chic;  
    pf.r_zlb = ones(length(pf.r),1);

    display(char('Policy Functions Calculated'))

    save('pf_cheb.mat','pf','O','P','S','G','C');
end


%--------------------------------------------------------------------------
% Construct CF Policy Functions
%--------------------------------------------------------------------------

pf.cf.c = zeros(G.griddim);
pf.cf.inf = zeros(G.griddim);
pf.cf.pitilde = zeros(G.griddim);
pf.cf.n = zeros(G.griddim);
pf.cf.y = zeros(G.griddim);
pf.cf.w = zeros(G.griddim);
pf.cf.r = zeros(G.griddim);

for i = 1:G.nodes
    if pf.r(i) >= 1 
        pf.cf.c(i) = pf.c(i);
        pf.cf.inf(i) = pf.inf(i);
        pf.cf.pitilde(i) = pf.pitilde(i);
        pf.cf.n(i) = pf.n(i);
        pf.cf.y(i) = pf.y(i);
        pf.cf.w(i) = pf.w(i);
        pf.cf.r(i) = pf.r(i);
    else
        pf.cf.c(i) = pf.c_zlb(i);
        pf.cf.inf(i) = pf.inf_zlb(i);
        pf.cf.pitilde(i) = pf.pitilde_zlb(i);
        pf.cf.n(i) = pf.n_zlb(i);
        pf.cf.y(i) = pf.y_zlb(i);
        pf.cf.w(i) = pf.w_zlb(i);
        pf.cf.r(i) = pf.r_zlb(i);
    end
end 

%--------------------------------------------------------------------------
% Plot Policy Functions
%-------------------------------------------------------------------------- 
zlb = find(pf.r < 1,1); %cell where ZLB binds for increasing \delta levels

X = [(pf.cf.c-S.c)./S.c, 400*(pf.cf.inf-1), (pf.cf.y-S.y)./S.y, 400*(pf.cf.r-1)];
% disp(400*(pf.cf.r-1));
title_fig = {'Consumption (S)', 'Inflation Rate (S)', 'Output (S)', 'Interest Rate (S)'};
title_fig_d = {'Consumption (D)', 'Inflation Rate (D)', 'Output (D)', 'Interest Rate (D)'};

ylim = [-.2 .1;-15 12;-.05 .15;0 20];

fig(1) = figure(1);
for i = 1:length(title_fig)
    subplot(2,2,i)
    box on
    hold on
    grid on
    plot(G.del_grid,X(:,i),'k','LineWidth',2)
    xlabel('\delta','FontSize',16)
    title(title_fig{i},'FontSize',16)
    set(gca,'XLim',[1-4*bound 1+4*bound],'FontSize',16)
    line([G.del_grid(zlb,1) G.del_grid(zlb,1)],get(gca,'YLim'),'Color','b','LineStyle','--','LineWidth',1.5)
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','Pfs_ZLBEconomy.eps');
