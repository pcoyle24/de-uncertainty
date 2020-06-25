% Fixed point iteration (Chebyshev polynomial basis):
% Canonical New Keynesian Model (Rotemberg Pricing) with Capital 
%   -Endogenous MP rule responds to current inflation (deterministic)
%   -Endogenous tax rate rule responds to lagged debt (stochastic)
%
% Unless otherwise noted, this script and the functions it depends on are
% authored by:
%   Alexander Richter (arichter@auburn.edu)
%   Auburn University - Department of Economics
%   Nathaniel Throckmorton (nathrock@indiana.edu)
%   Indiana University - Department of Economics
%
% Copyright 2010-2013: You may freely copy and alter the following code,
% (including any scripts, functions and documentation authored by us that
% this script depends on) under the following conditions:
%    1) This code and any modifications to it are not sold or
%       traded in exchange for any goods or services
%    2) Credit is given to the original authors upon redistribution of
%       the original or modified code
%    3) This copyright notice is included on the original code and
%       subsequent modifications


clear all
close all
clc

tstart = tic;                           % Job timer start
%--------------------------------------------------------------------------
% Initialize Policy Functions
%--------------------------------------------------------------------------
% Load parameters, steady state and grids
P = parameters;
S = steadystate(P);

% Specify grid options
O.abound = [0.95*S.a 1.05*S.a];
O.bbound = [0.95*S.b 1.05*S.b];
O.kbound = [0.95*S.k 1.05*S.k];
O.fpbound = [-0.05 0.05];
O.a_pts = 7;
O.b_pts = 7;
O.k_pts = 7;
O.fp_pts = 7;
O.e_pts = 10;

% Specify Chebyshev polynomial orders (must <= pts in grid)
O.n1 = 4;
O.n2 = 4;
O.n3 = 4;
O.n4 = 4;

% Load discretized state space
G = grids_cheb(O,P,S);

% Initialize Chebyshev polynomial parameters
C = chebpoly(G.nl.griddim);

% Specify initial conjecture
%  guess: linear policies
%  current: current nonlinear policy functions
O.loadpf = 'guess';

% Retrieve initial policy functions
pf = guess(O,P,S,G); 

% pf1 = guess_flat(O,P,S,G);
%--------------------------------------------------------------------------
% Initialize algorithm
%--------------------------------------------------------------------------
% Calculate initial Chebychev coefficients
[C.An C.Api C.Ak] = chebweights43(  O.n1,O.n2,O.n3,O.n4, ...
                                       pf.n,pf.pi,pf.k, ...
                                       C.T,C.P,C.X);
% Vectorize grids and policy function
aMat = G.nl.a_gr(:,:,:,:,ones(1,1,O.e_pts));
bMat = G.nl.b_gr(:,:,:,:,ones(1,1,O.e_pts));
kMat = G.nl.k_gr(:,:,:,:,ones(1,1,O.e_pts));
fpMat = G.nl.fp_gr(:,:,:,:,ones(1,1,O.e_pts));
eMat = G.e_nodes(:,ones(O.a_pts,1),ones(O.b_pts,1),ones(O.k_pts,1),ones(O.fp_pts,1));
eMat = permute(eMat,[2 3 4 5 1]);
eWeight = G.e_weight(:,ones(O.a_pts,1),ones(O.b_pts,1),ones(O.k_pts,1),ones(O.fp_pts,1));
eWeight = permute(eWeight,[2 3 4 5 1]);
nodes = numel(kMat);
%--------------------------------------------------------------------------
% Fixed point iteration/Chebyshev polynomial basis algorithm
%--------------------------------------------------------------------------
it = 1;                                 % Iteration Counter
converged = 0;                          % Convergence Flag
while converged == 0        
    istart = tic;                       % Iteration timer start    
    
    % Map policies over the shocks    
    nMat = pf.n(:,:,:,:,ones(1,O.e_pts));
    piMat = pf.pi(:,:,:,:,ones(1,O.e_pts));
    kpMat = pf.k(:,:,:,:,ones(1,O.e_pts));
    
    %----------------------------------------------------------------------
    % Solve for current period variables
    %---------------------------------------------------------------------- 
    % Interest rate rule
    rMat = S.r*(piMat/P.pi).^P.phi;
    % Tax rule
    tauMat = P.tau*(bMat/S.b).^P.gam.*exp(fpMat);
    % Production function    
    yMat = kMat.^P.alpha.*nMat.^(1-P.alpha);
    % Investment 
    iMat = kpMat - (1-P.delta)*kMat;         
    % Aggregate resource constraint    
    cMat = yMat.*(1-(P.varphi*(piMat/P.pi-1).^2)/2) - iMat - S.g;
    % FOC Labor
    wMat = S.chi*nMat.^P.eta.*cMat.^P.sigma./(1-tauMat);    
    % Combine: Firm FOCs
    rkMat = wMat.*nMat*P.alpha./(kMat*(1-P.alpha));    
    % FOC Money
    mMat = S.xi^(1/P.kappa)*((rMat-1)./rMat).^(-1/P.kappa).*cMat.^(P.sigma/P.kappa);        
    % Government budget constraint
    bpMat = aMat./piMat + S.g - tauMat.*(wMat.*nMat + rkMat.*kMat) - mMat;           
    % Marginal cost definition
    psiMat = wMat.^(1-P.alpha).*rkMat.^P.alpha/((1-P.alpha)^(1-P.alpha)*P.alpha^P.alpha);        
    % Government liabilities definition
    apMat = mMat + rMat.*bpMat;           
    %----------------------------------------------------------------------
    % Calculate updated policy function (t-1 coefficients) 
    %----------------------------------------------------------------------
    [npMat pipMat kppMat] = ...
    allcheb431_F(  O.abound,O.bbound,O.kbound,O.fpbound, ...
                    apMat,bpMat,kpMat,eMat, ...
                    C.An,C.Api,C.Ak,C.T,C.P);
    %----------------------------------------------------------------------
    % Solve for next period variables
    %----------------------------------------------------------------------
    % Tax rate
    taupMat = P.tau*(bpMat/S.b).^P.gam.*exp(eMat);
    % Production function    
    ypMat = kpMat.^P.alpha.*npMat.^(1-P.alpha);
    % Investment
    ipMat = kppMat-(1-P.delta)*kpMat;
    % Aggregate resource constraint 
    cpMat = ypMat.*(1-(P.varphi*(pipMat./P.pi-1).^2)/2)-S.g-ipMat;
    % FOC Labor
    wpMat = S.chi*(npMat.^P.eta).*(cpMat.^P.sigma)./(1-taupMat);
    % Combine: Firm FOCs
    rkpMat = wpMat.*npMat.*P.alpha./(kpMat*(1-P.alpha));
    % Stochastic discount factor
    qMat = P.beta*(cMat./cpMat).^P.sigma;        
    %----------------------------------------------------------------------
    % Evaluate expectations (GH Quadrature) 
    %----------------------------------------------------------------------
    % Consumption Euler equation
    MUcpMat = cpMat.^(-P.sigma);
    Econs =  pi^(-.5)*eWeight.*(P.beta.*MUcpMat.*((1-taupMat).*rkpMat+1-P.delta));
    Econs = sum(Econs,5);
    % Firm pricing equation
    Efp = pi^(-.5)*eWeight.*(qMat.*(pipMat./P.pi-1).*(ypMat./yMat).*(pipMat./P.pi));
    Efp = sum(Efp,5);
    % Bond Euler Equation
    Ebond = pi^(-.5)*eWeight.*(qMat./pipMat);
    Ebond = sum(Ebond,5);   
    %----------------------------------------------------------------------
    % Solve for policy functions
    %----------------------------------------------------------------------
    % Bond Euler equation
    pf_r = 1./Ebond; %exp
    % Interest rate rule
    pf.pi = P.pi*(pf_r/S.r).^(1/P.phi);
    % Firm Pricing
    pf_psi = (P.varphi*(pf.pi/P.pi-1).*pf.pi/P.pi-(1-P.theta)-P.varphi*Efp)/P.theta;
    % Consumption euler
    pf_c = Econs.^(-1/P.sigma);  
    % Combine: FOC labor, firm FOC labor, and production function
    pf.n = ((1-tauMat(:,:,:,:,1)).*(1-P.alpha).*G.nl.k_gr.^P.alpha.*pf_psi./(S.chi*pf_c.^P.sigma)).^(1/(P.eta+P.alpha));        
    % Production function
    pf_y = G.nl.k_gr.^P.alpha.*pf.n.^(1-P.alpha);
    % Aggregate resource constraint
    pf_i = pf_y.*(1-(P.varphi*(pf.pi/P.pi-1).^2)/2) - pf_c - S.g;
    % Investment
    pf.k = pf_i + (1-P.delta)*G.nl.k_gr;                        
     
    % New policy function coefficients
    [An_up Api_up Ak_up] = chebweights43(O.n1,O.n2,O.n3,O.n4, ...
                                            pf.n,pf.pi,pf.k, ...
                                            C.T,C.P,C.X);

    % Coefficient distances
    dist_An = abs(An_up - C.An);
    dist_Api = abs(Api_up - C.Api);
    dist_Ak = abs(Ak_up - C.Ak);
    
    % Maximum distance
    dist_max = max([dist_An(:)' dist_Api(:)' dist_Ak(:)']);
    
    % Update policy function coefficients    
    C.An = P.lambda*An_up + (1-P.lambda)*C.An;
    C.Api = P.lambda*Api_up + (1-P.lambda)*C.Api;
    C.Ak = P.lambda*Ak_up + (1-P.lambda)*C.Ak;
    
    % Calculate updated policy function                  
    [pf.n pf.pi pf.k] = ...
    allcheb431_F(  O.abound,O.bbound,O.kbound,O.fpbound, ...
                    aMat,bMat,kMat,fpMat, ...
                    C.An,C.Api,C.Ak,C.T,C.P);   

    % Resize policy functions    
    pf.n = pf.n(:,:,:,:,1);
    pf.pi = pf.pi(:,:,:,:,1);
    pf.k = pf.k(:,:,:,:,1);

    % Save policy functions
    save('pf.mat','pf','O','P','S','G','C')
    
    % Iteration Information
    it = itinfo(istart,tstart,1,it,dist_max);

    % Check convergence criterion
    if dist_max < P.tol; converged = 1; end;
end