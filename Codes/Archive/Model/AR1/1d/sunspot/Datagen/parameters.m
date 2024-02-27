function P = parameters

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
P.alpha    = 0;
P.iota     = 1;
P.pi_targ  = 2/400 + 1;

% Monetary Policy 
P.phi_pi     = 2;       % Inflation coefficient: active interest rate rule
P.phi_y      = 0;       % Inflation coefficient: active interest rate rule

% Fiscal Policy
P.rho     = 0.8;       % Debt coefficient: passive tax rule
P.sigma_grid   = [0.15/100, 0.25/100, 0.27325/100, 0.32/100, 0.335/100];


% Transition Probabilities 
P.Ps_grid      = [0.995,1]; % Probability of transitioning from SSS to SSS 
P.Pd_grid      = [1, 0.975]; % Probability of transitioning from DSS to DSS
