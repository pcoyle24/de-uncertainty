function P = parameters

% P = parameters
%   Defines model parameters
% Output:
%   P : Structure of model parameters

% Household (Annual Calibration)
P.beta     = 0.9615;   % Discount factor
P.sigma    = 1;        % Constant of relative risk aversion
P.delta    = 0.1;      % Depreciation rate of capital
P.eta      = 1;        % Inverse frish elasticity of labor supply
P.kappa    = 1;        % Inverse interest (semi) elasticity of money demand
P.theta    = 7.6667;   % Price elasticity of demand
P.vel      = 3.8;      % Money velocity
P.n        = 0.33;     % Steady state labor

% Firm
P.alpha    = 0.33;     % Cost share of capital
P.varphi   = 10;       % Rotemberg adjustment cost coefficient

% Monetary Policy 
P.pi      = 1.02;      % Inflation target
P.phi     = 1.5;       % Inflation coefficient: active interest rate rule

% Fiscal Policy
P.gy      = 0.17;      % Government spending share of GDP
P.tau     = 0.21;      % Steady state tax rate
P.gam     = 0.2;       % Debt coefficient: passive tax rate rule
P.sig_fp  = 0.001;     % Standard deviation of tax shock 

% Algorithm
P.tol     = 1e-10;     % Convergence criterion
P.lambda  = 0.5;       % Weight on updated coeffiecients 