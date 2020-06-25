% Projection Method (Chebyshev basis):
% RBC Model with Stochastic AR(1) Productivity Processes
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

clear;
clc;
tstart = tic;                           % Job timer start
%--------------------------------------------------------------------------
% Initialize Policy Functions
%--------------------------------------------------------------------------
% Load parameters, steady state and grids
P = parameters;
S = steadystate(P);

% Specify grid options
O.kbound = [S.k*0.97 S.k*1.03];
O.zbound = [0.97 1.03];
O.k_pts = 21;
O.z_pts = 7;
O.e_pts = 10;

% Chebyshev polynomial order in each dimension (must <= pts in grid)
O.n1 = 5;
O.n2 = 5;

% Load discretized state space
G = grids_cheb(O,P);

% Initialize Chebyshev polynomial parameters
C = chebpoly(G.griddim);

% Specify initial conjecture
%  guess: linear policies
%  current: current nonlinear policy functions
O.loadpf = 'guess';
O.db = dbstack;

% Retrieve initial policy functions
pf = guess(O,P,S,G);
% pf.n = ones(G.griddim)*(P.n);
%--------------------------------------------------------------------------
% Initialize algorithm
%--------------------------------------------------------------------------
% Calculate initial Chebyshev coefficients          
C.An = chebweights21(O.n1,O.n2,pf.n,C.T,C.P,C.X);

% Reshape state and shocks
kMat = G.k_gr(:,:,ones(1,1,O.e_pts));
zMat = G.z_gr(:,:,ones(1,1,O.e_pts));
eMat = G.e_nodes(:,ones(O.k_pts,1));
eMat = eMat(:,:,ones(O.z_pts,1));
eMat = permute(eMat,[2 3 1]);
eWeight = G.e_weight(:,ones(O.k_pts,1));
eWeight = eWeight(:,:,ones(O.z_pts,1));
eWeight = permute(eWeight,[2 3 1]);
nodes = numel(kMat);
griddim = numel(kMat);
%--------------------------------------------------------------------------
% Least squares projection algorithm
%--------------------------------------------------------------------------
it = 1;                                 % Iteration Counter
converged = 0;                          % Convergence Flag
while converged == 0
    istart = tic;                       % Iteration timer start    
    
    % Replicate policy function over shocks
    nMat = pf.n(:,:,ones(1,O.e_pts));    
    %----------------------------------------------------------------------
    % Solve for current period variables
    %----------------------------------------------------------------------    
    yMat = zMat.*kMat.^P.alpha.*nMat.^(1-P.alpha);
    cMat = ((1-P.alpha)*yMat./(S.chi*nMat.^(1+P.eta))).^(1/P.sigma); 
    iMat = yMat - cMat;    
    kpMat = iMat + (1-P.delta)*kMat; 
    wMat = S.chi*cMat.^P.sigma.*nMat.^P.eta;      
    zpMat = (1-P.rho)*P.zbar + P.rho*zMat + eMat;
    %----------------------------------------------------------------------
    % Calculate updated policy function (t-1 coefficients) 
    %----------------------------------------------------------------------
    npMat = allcheb211_F(O.kbound,O.zbound,kpMat,zpMat,C.An,C.T,C.P);    
    %----------------------------------------------------------------------
    % Solve for next period variables
    %----------------------------------------------------------------------
    ypMat = zpMat.*kpMat.^P.alpha.*npMat.^(1-P.alpha);
    cpMat = ((1-P.alpha)*ypMat./(S.chi*npMat.^(1+P.eta))).^(1/P.sigma); 
    rkpMat = P.alpha*ypMat./kpMat;
    %----------------------------------------------------------------------
    % Evaluate expectations (GH Quadrature) 
    %----------------------------------------------------------------------
    MUcpMat = cpMat.^(-P.sigma);
    Econs =  pi^(-.5)*eWeight.*(P.beta.*MUcpMat.*(rkpMat+1-P.delta));
    Econs = sum(Econs,3);
    %----------------------------------------------------------------------
    % Solve for policy functions
    %----------------------------------------------------------------------
    % Consumption
     pf_c = Econs.^(-1/P.sigma);
    % Labor
    pf.n = (wMat(:,:,1)./(S.chi*pf_c.^P.sigma)).^(1/P.eta);
    
    % New policy function coefficients        
    An_up = chebweights21(O.n1,O.n2,pf.n,C.T,C.P,C.X);  
    
    % Coefficient distances   
    dist_An = abs(An_up - C.An);
    
    % Maximum distance
    dist_max = max(dist_An(:));
    
    % Update policy function coefficients    
    C.An = P.lambda*An_up + (1-P.lambda)*C.An;
    
    % Calculate updated policy function
    pf.n = allcheb211_F(O.kbound,O.zbound,kMat,zMat,C.An,C.T,C.P);
    
    % Save policy functions
    save('pf_flat.mat','pf','O','P','S','G','C')

    % Iteration Information
    it = itinfo(istart,tstart,1,it,dist_max);
    
    % Check convergence criterion
    if dist_max < P.tol; converged = 1; end;
end