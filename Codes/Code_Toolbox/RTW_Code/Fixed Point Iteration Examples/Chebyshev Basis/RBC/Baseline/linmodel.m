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
G0(j,V.n) = -P.eta;
G0(j,V.c) = -P.sigma;
%---------------------------------------------
j = j+1;% [Consumption Euler Equation]
%---------------------------------------------
G0(j,V.c) = P.sigma;
G0(j,V.rk) = -(1-P.beta*(1-P.delta));
G1(j,V.c) = P.sigma;
Pi(j,V.fec) = P.sigma; 
Pi(j,V.ferk) = -(1-P.beta*(1-P.delta));
%---------------------------------------------
j = j+1;% [Aggregate Resource Constraint]
%---------------------------------------------
G0(j,V.c) = S.c;
G0(j,V.y) = -S.y;
G0(j,V.i) = S.i;
%---------------------------------------------
j = j+1;% [Production Function] 
%---------------------------------------------
G0(j,V.y) = 1;
G0(j,V.z) = -1;
G0(j,V.n) = -(1-P.alpha);
G1(j,V.k) = P.alpha;
%---------------------------------------------
j = j+1;% [Law of Motion for Capital]
%---------------------------------------------
G0(j,V.i) = P.delta;
G0(j,V.k) = -1;
G1(j,V.k) = -(1-P.delta);
%---------------------------------------------
j = j+1;% [Firm FOC Labor]
%---------------------------------------------
G0(j,V.w) = 1;
G0(j,V.y) = -1;
G0(j,V.n) = 1;
%---------------------------------------------
j = j+1;% [Firm FOC Capital]
%---------------------------------------------
G0(j,V.rk) = 1;
G0(j,V.y) = -1;
G1(j,V.k) = -1;
%---------------------------------------------
j = j+1;% [Productivity]
%---------------------------------------------
G0(j,V.z) = 1;
G1(j,V.z) = P.rho;
Psi(j,V.epsz) = 1;

%---------------------------------------------
%   Solve Linear Model
%---------------------------------------------
[T,~,M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);           