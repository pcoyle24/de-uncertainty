function P = parameters

% Parameters: Specifies model parameters
%   P = parameters creates a structure of model parameters

% Household (Quarterly Calibration)
P.beta    = 0.99;     % Discount factor
P.delta   = 0.025;    % Depreciation rate of capital
P.n       = 0.33;  	  % Steady state labor
P.sigma   = 2;        % Risk aversion
P.eta     = 2;        % Inverse Frisch elasticity of labor supply

% Firm
P.alpha   = 0.3;      % Cost share of capital

% Productivity
P.zbar    = 1;        % Productivity process intercept
P.rho     = 0.9;      % AR coefficient: productivity process  
P.sig_eps = 0.001;    % Productivity Shock Std. Dev.

% Algorithm
P.tol     = 1e-10;    % Convergence criterion
P.lambda  = 0.50;     % Weight on updated coeffiecients 