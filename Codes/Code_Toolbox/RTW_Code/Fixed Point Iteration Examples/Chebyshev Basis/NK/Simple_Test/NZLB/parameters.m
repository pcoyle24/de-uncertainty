function P = parameters

% P = parameters
%   Defines model parameters
% Output:
%   P : Structure of model parameters

% Household (Annual Calibration)
P.beta     = (1/1.0075);   % Discount factor
P.chic     = 1;            % Constant of relative risk aversion
P.chin     = 1;            % Constant of relative risk aversion
P.theta    = 10;           % CES technology production


% Firm
P.varphi   = 175;       % Rotemberg adjustment cost coefficient

% Monetary Policy 
P.phi_pi     = 1.5;       % Inflation coefficient: active interest rate rule

% Fiscal Policy
P.rho     = 0.8;       % Debt coefficient: passive tax rule
P.sigma   = 0/100;   % Standard deviation of demand shock 
% P.sigma   =0;   % Standard deviation of demand shock 
% Algorithm
P.tol     = 1e-10;     % Convergence criterion
P.lambda  = 1;     % Weight on updated coeffiecients 