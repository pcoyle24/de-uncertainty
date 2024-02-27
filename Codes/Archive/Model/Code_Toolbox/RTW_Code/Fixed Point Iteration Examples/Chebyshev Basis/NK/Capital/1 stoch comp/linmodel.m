function [T, M, eu] = linmodel(P,S,V)

% [T, M, eu] = linmodel(P,S,V)
%   Solves the log-linear model with GENSYS
% Inputs:
%     P     :   Structure of parameters
%     S     :   Structure of steady state values
%     V     :   Structure of variable locations and names
% Output:
%     T     :   Transition matrix
%     M     :   Impact matrix
%     eu    :   [existence, uniqueness]

%---------------------------------------------
%   Initialize GENSYS components
%---------------------------------------------
G0 = zeros(V.nvar);
G1 = zeros(V.nvar);
Psi = zeros(V.nvar,V.nshock);
Pi = zeros(V.nvar,V.nfore);
CC = zeros(V.nvar,1);
j = 0;
%---------------------------------------------
j=j+1;%	[FOC Labor]
%---------------------------------------------
G0(j,V.w) = 1;
G0(j,V.c) = -P.sigma;
G0(j,V.n) = -P.eta;
G0(j,V.tau) = -P.tau/(1-P.tau);
%---------------------------------------------
j = j+1;% [Consumption Euler Equation]
%---------------------------------------------
G0(j,V.c) = P.sigma;
G0(j,V.rk) = -(1-P.beta*(1-P.delta));
G0(j,V.tau)= P.tau*(1-P.beta*(1-P.delta))/(1-P.tau);
G1(j,V.c) = P.sigma;
Pi(j,V.fec) = P.sigma; 
Pi(j,V.ferk) = -(1-P.beta*(1-P.delta));
Pi(j,V.fetau)= P.tau*(1-P.beta*(1-P.delta))/(1-P.tau); 
%---------------------------------------------
j = j+1;% [FOC Money]
%---------------------------------------------
G0(j,V.m) = 1;
G0(j,V.c) = -P.sigma/P.kappa;
G0(j,V.r) = 1/((S.r-1)*P.kappa);
%---------------------------------------------
j = j+1;% [FOC Bond]
%---------------------------------------------
G0(j,V.c) = P.sigma;
G0(j,V.pi) = 1;
G1(j,V.r) = 1;
G1(j,V.c) = P.sigma;
Pi(j,V.fec) = P.sigma; 
Pi(j,V.fepi) = 1; 
%---------------------------------------------
j = j+1;% [ARC]
%---------------------------------------------
G0(j,V.c) = S.c;
G0(j,V.y) = -S.y;
G0(j,V.i) = S.i;
%---------------------------------------------
j = j+1;% [Production Function] 
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.n) = -(1-P.alpha);
G1(j,V.k) = P.alpha;
%---------------------------------------------
j = j+1;% [Firm Pricing] 
%---------------------------------------------
G0(j,V.pi) = P.beta;
G1(j,V.pi) = 1;
G1(j,V.psi) = (1-P.theta)/P.varphi;
Pi(j,V.fepi) = P.beta; 
%---------------------------------------------
j = j+1;% [Marginal Cost Definition]
%---------------------------------------------
G0(j,V.psi) = 1;
G0(j,V.w) = -(1-P.alpha);
G0(j,V.rk) = -P.alpha;
%---------------------------------------------
j = j+1;% [GBC] 
%---------------------------------------------
G0(j,V.tr) = S.tr;
G0(j,V.m) = S.m;
G0(j,V.b) = S.b; 
G0(j,V.pi) = S.a/P.pi; 
G0(j,V.a) = -S.a/P.pi;
%---------------------------------------------
j = j+1;% [Tax Process]
%---------------------------------------------
G0(j,V.tau) = 1;
G1(j,V.b) = P.gam;
Psi(j,V.epsfp) = 1;
%---------------------------------------------
j = j+1;% [Interest Rate Rule]
%---------------------------------------------
G0(j,V.r) = 1; 
G0(j,V.pi) = -P.phi;
%---------------------------------------------
j = j+1;% [Law of Motion for Capital]
%---------------------------------------------
G0(j,V.i) = P.delta;
G0(j,V.k) = -1;
G1(j,V.k) = -(1-P.delta);
%---------------------------------------------
j = j+1;% [Firm Pricing] 
%---------------------------------------------
G0(j,V.w) = 1; 
G0(j,V.n) = 1;
G0(j,V.rk) = -1;
G1(j,V.k) = 1;
%---------------------------------------------
j = j+1;% [Government Liabilities Definition] 
%---------------------------------------------
G0(j,V.a) = S.a; 
G1(j,V.m) = S.m;
G1(j,V.r) = S.r*S.b;
G1(j,V.b) = S.r*S.b;
%---------------------------------------------
j = j+1;% [Tax Revenues] 
%---------------------------------------------
G0(j,V.tr) = S.tr;
G0(j,V.tau) = -S.tr; 
G0(j,V.w) = -P.tau*S.w*P.n; 
G0(j,V.n) = -P.tau*S.w*P.n;
G0(j,V.rk) = -P.tau*S.rk*S.k;
G1(j,V.k) = P.tau*S.rk*S.k;

%---------------------------------------------
%   Solve Linear Model
%---------------------------------------------
[T,~,M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);        