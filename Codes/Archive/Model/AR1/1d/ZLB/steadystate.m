function S = steadystate(P,i)
% S = steadystate(P)
%   Computes the deterministic steady state with a nonlinear solver
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

S.pi_targ = P.pi_targ(i);

x0 = [0.9 0.9 0.9 -67];
options = optimset('display','none','MaxFunEvals',100000,'MaxIter',1000,'TolFun',1e-32,'TolX',1e-32,'FunValCheck','On');

%% Standard Steady State
func = @(x) steady(x,P,i);
% func = @(x) steady_deflationaryE(x,P,i);
out = fsolve(func,x0,options);

S.c    = out(1);
S.inf  = out(2);
S.n    = out(3);
S.v    = out(4);
S.y    = S.n;
S.w    = S.n^P.chin*S.c^P.chic;
S.r    = S.pi_targ/P.beta*((S.inf/S.pi_targ)^(P.phi_pi)*(S.y/S.y)^(P.phi_y));

S.del    = 1;
S.r_zlb  = 1;


%% Demand Shock Steady State
S.del = 1; 

