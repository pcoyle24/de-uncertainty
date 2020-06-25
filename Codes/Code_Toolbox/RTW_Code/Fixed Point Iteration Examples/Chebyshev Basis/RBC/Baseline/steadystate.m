function S = steadystate(P)

% S = steadystate(P)
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% Consumption Euler Equation
S.rk = (1/P.beta + P.delta - 1);
% Firm FOC Capital
ky = P.alpha/S.rk;
% Prod. Func.
yn = ky^(P.alpha/(1-P.alpha));
% Ratios
S.y = yn*P.n;
S.k = ky*S.y;
% Investment
S.i = P.delta*S.k;
% Aggregate resouce constraint
S.c = S.y - S.i;
% Firm FOC Labor
S.w = (1-P.alpha)*S.y/P.n;
% Household FOC labor
S.chi = S.w/(S.c^P.sigma*P.n^P.eta);