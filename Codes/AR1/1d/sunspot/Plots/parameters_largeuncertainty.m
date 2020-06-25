function P = parameters_largeuncertainty

% P = parameters
%   Defines model parameters
% Output:
%   P : Structure of model parameters

% Household (Annual Calibration)
P.beta     = (1/1.0025);   % Discount factor
P.chic     = 1;            % Constant of relative risk aversion
P.chin     = 1;            % Constant of relative risk aversion
P.theta    = 11;           % CES technology production


% Firm
P.varphi   = 500;       % Rotemberg adjustment cost coefficient
P.tau      = 1/P.theta; 
P.alpha    = 1;
P.iota     = 1;
P.pi_targ  = 2/400 + 1;

% Monetary Policy 
P.phi_pi     = 16;       % Inflation coefficient: active interest rate rule
P.phi_y      = 0;       % Inflation coefficient: active interest rate rule

% Fiscal Policy
P.rho     = 0.8;       % Debt coefficient: passive tax rule
% P.sigma   = 0.2745/100;   % Standard deviation of demand shock
P.sigma   = 0.419225/100;%0.493/100;   % Standard deviation of demand shock
P.bound = (P.sigma^2/(1-P.rho^2))^0.5;



% Transition Probabilities 
P.Ps      = 1; % Probability of transitioning from SSS to SSS 
P.Pd      = 1; % Probability of transitioning from DSS to DSS
