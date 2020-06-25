%--------------------------------------------------------------------------
%File Name: script.m
%Author: Philip Coyle
%Date Created: 02/06/2018
%cd P:\RA_Work\Taisuke_Nakata\Zero_Lower_Bound\DeflationaryRegime\Cheb_play\Code_Toolbox\RTW_Code\Fixed Point Iteration Examples\Chebyshev Basis\NK\Simple_Test\NZLB
%--------------------------------------------------------------------------

% Projection Method (Chebyshev basis):
% Canonical New Keynesian Model (Rotemberg Pricing) with AR1 Preference
% Shock


clear all;
close all;
clc;
tstart = tic;                           % Job timer start

%--------------------------------------------------------------------------
% Initialize Policy Functions
%--------------------------------------------------------------------------
% Load parameters, steady state and grids
P = parameters;
S = steadystate(P);

% Specify grid options
O.delbound = [0.975*S.del 1.025*S.del];
O.del_pts = 11;
O.e_pts = 10;

% Chebyshev polynomial order in each dimension (must <= pts in grid)
O.n1 = 5;


% Load discretized state space
G = grids_cheb(O,P);

% Initialize Chebyshev polynomial parameters
C = chebpoly(G,O);

% % Retrieve initial policy functions
% pf = guess_flat(O,P,S,G); 

% Specify initial conjecture
O.loadpf = 'guess';

% Retrieve initial policy functions
pf = guess1(O,P,S,G); 


%--------------------------------------------------------------------------
% Initialize algorithm
%--------------------------------------------------------------------------
% Calculate initial Chebyshev coefficients          
C.Ac = chebweights11(O.n1,pf.c,C.T,C.P,C.X);
% pf_c = testchebbuild(O.delbound,C.Ac,C.T,C.P,O.e_pts,G.del_gr);
C.Ainf = chebweights11(O.n1,pf.inf,C.T,C.P,C.X);
C.Ar = chebweights11(O.n1,pf.r,C.T,C.P,C.X);

% Reshape state and shockse_pts
del_today = G.del_gr(:,ones(1,1,O.e_pts));
% del_today = ones(size(G.del_gr(:,ones(1,1,O.e_pts))));

del_shock = G.e_nodes(:,ones(O.del_pts,1));
del_shock = permute(del_shock,[2 1]);

e_weight = G.e_weight(:,ones(O.del_pts,1));
e_weight = permute(e_weight,[2 1]);


it = 1;                                 % Iteration Counter
converged = 0;                          % Convergence Flag
while converged == 0
    istart = tic;                       % Iteration timer start    
    
    pf.n = pf.c./(1-P.varphi/2*(pf.inf - 1).^2);
    % Replicate policy function over shocks
    c_today = pf.c(:,ones(1,O.e_pts));    
    inf_today = pf.inf(:,ones(1,O.e_pts));  
    %----------------------------------------------------------------------
    % Solve for current period variables
    %----------------------------------------------------------------------    
    n_today = c_today./(1-P.varphi/2*(inf_today - 1).^2);
    r_today = 1/P.beta*inf_today.^P.phi_pi;
    
    del_tomorrow = P.rho*(del_today - 1) + 1 + del_shock;
    %----------------------------------------------------------------------
    % Calculate updated policy function (t-1 coefficients) 
    %----------------------------------------------------------------------
    c_tomorrow = allcheb111_F(O.delbound,del_tomorrow,C.Ac,C.T,C.P);
    inf_tomorrow = allcheb111_F(O.delbound,del_tomorrow,C.Ainf,C.T,C.P);
    
    %----------------------------------------------------------------------
    % Solve for next period variables
    %----------------------------------------------------------------------
    n_tomorrow = c_tomorrow./(1-P.varphi/2*(inf_tomorrow - 1).^2);
    %----------------------------------------------------------------------
    % Evaluate expectations (GH Quadrature) 
    %----------------------------------------------------------------------

    exp_ee_c =  P.beta.*r_today.*pi^(-.5).*e_weight.*(del_today.*c_tomorrow.^(-P.chic).*inf_tomorrow.^(-1));
    exp_ee_c = sum(exp_ee_c,2);
    
%     exp_ee_r =  (pi^(-.5)*e_weight.*(P.beta.*del_today.*(c_today./c_tomorrow).^(P.chic).*inf_tomorrow.^(-1)));
%     exp_ee_r = sum(exp_ee_r,2);
%     
    
    exp_pc = pi^(-.5)*e_weight.*(P.beta.*del_today.*((n_tomorrow./n_today).*(c_today./c_tomorrow).^(P.chic)).*P.varphi.*(inf_tomorrow-1).*inf_tomorrow);
    exp_pc = sum(exp_pc,2);
%     exp_pc = exp_pc(:,ones(1,1,O.e_pts)); 
    
    pc_rhs = exp_pc + (1-P.theta) + P.theta.*(pf.n.^(P.chin).*pf.c.^(P.chic));
    
    a = ones(length(pc_rhs),1)*P.varphi;
    b = -a;
    pi_root = quadroot(a,b,-pc_rhs);
    


%     %Back out via pf_n
%     exp_pc = pi^(-.5)*e_weight.*(P.beta.*del_today.*((n_tomorrow./n_today).*(c_today./c_tomorrow).^(P.chic)).*P.varphi.*(inf_tomorrow-1).*inf_tomorrow);
%     exp_pc = sum(exp_pc,2);
%     exp_pc = exp_pc(:,ones(1,1,O.e_pts));
%     pc_rhs_terms = P.varphi.*(inf_today - 1).*inf_today -(1-P.theta);
% 
% 
%     LHS = c_today.^(-P.chic).*(pc_rhs_terms - exp_pc);
%     
%     pf_n = (LHS./(P.theta)).^(1/P.chin);
    
    
    %----------------------------------------------------------------------
    % Solve for policy functions (Back them out)
    %----------------------------------------------------------------------
    % Consumption
    pf.c = (1./(exp_ee_c)).^(1/P.chic);
    
    % Inflation 
    pf.inf = max(pi_root,[],2);
    
    % Interst Rate
    pf.r = 1/P.beta.*(pf.inf).^P.phi_pi;
    
    % New policy function coefficients        
    Ac_up = chebweights11(O.n1,pf.c,C.T,C.P,C.X);
    Ainf_up = chebweights11(O.n1,pf.inf,C.T,C.P,C.X);
    Ar_up = chebweights11(O.n1,pf.r,C.T,C.P,C.X);
    
    % Coefficient distances   
    diff_Ac = abs(Ac_up - C.Ac);
    diff_Ainf = abs(Ainf_up - C.Ainf);
    diff_Ar = abs(Ar_up - C.Ar);
    
    % Maximum distance
    dist_max = max([diff_Ac(:)', diff_Ainf(:)', diff_Ar(:)']);
    
    % Update policy function coefficients    
    C.Ac = P.lambda*Ac_up + (1-P.lambda)*C.Ac;
    C.Ainf = P.lambda*Ainf_up + (1-P.lambda)*C.Ainf;
    C.Ar = P.lambda*Ar_up + (1-P.lambda)*C.Ar;
    
%     % Calculate updated policy function
%     pf.c = allcheb111_F(O.delbound,del_today,C.Ac,C.T,C.P);
%     pf.inf = allcheb111_F(O.delbound,del_today,C.Ainf,C.T,C.P);
    
%     pf.c(:,1);
    
%     pf.r = 1/P.beta*(pf.inf.^(P.phi_pi));
    
    % Save policy functions
%     save('pf.mat','pf','O','P','S','G','C')

    % Iteration Information
    it = itinfo(istart,tstart,1,it,dist_max);
    
    % Check convergence criterion
    if dist_max < P.tol; converged = 1; end;
end