function S = steadystate(P)

% S = steadystate(P)
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% FOC bonds 
S.r = P.pi/P.beta;
% FOC capital
S.rk = (1/P.beta + P.delta - 1)/(1-P.tau);
% Firm pricing
S.psi = (P.theta-1)/P.theta;
% Marginal cost definition
S.w = (S.psi*(1-P.alpha)^(1-P.alpha)*P.alpha^P.alpha/S.rk^P.alpha)^(1/(1-P.alpha));
% Consolidated FOC firm
S.k = S.w*P.n*P.alpha/(S.rk*(1-P.alpha));
% Investment definition
S.i = P.delta*S.k;
% Production function
S.y = S.k^P.alpha*P.n^(1-P.alpha);
% Government spending (level)
S.g = P.gy*S.y;
% Aggregate resouce constraint
S.c = S.y - S.i - S.g;
% FOC labor
S.chi = S.w*(1-P.tau)/(P.n^P.eta*S.c^P.sigma);
% Velocity definition
S.m = S.c/P.vel;
% FOC Money
S.xi = S.m^P.kappa*(S.r-1)/(S.c^P.sigma*S.r);
% Tax revenue definition
S.tr = P.tau*S.w*P.n + P.tau*S.rk*S.k;
% Government budget constraint
S.b = (S.g - S.tr - S.m*(1-1/P.pi))/(1-S.r/P.pi);
% Government liabilities
S.a = S.m + S.r*S.b;