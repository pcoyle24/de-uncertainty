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
CC = zeros(V.nvar,1);
Psi = zeros(V.nvar,V.nshock);
Pi = zeros(V.nvar,V.nfore);
j = 0;
%---------------------------------------------
j=j+1;%	[Consumption Euler Equation]
%---------------------------------------------
G0(j,V.c)     = -P.chic;
G0(j,V.inf)   = -1;
G1(j,V.c)     = -P.chic;
G1(j,V.del)   = -1;
G1(j,V.r)     = -1;
Pi(j,V.feinf) = -1;
Pi(j,V.fec)   = -P.chic;
%---------------------------------------------
j = j+1;% [FOC Labor]
%---------------------------------------------
G0(j,V.inf)   = P.beta*P.varphi;
G1(j,V.c)     = P.chic*(1-P.theta);
G1(j,V.n)     = -(1-P.theta+P.theta*S.n^(P.chin)*S.c^(P.chic)*(1+P.chin));
G1(j,V.inf)   = P.varphi;
Pi(j,V.feinf) = P.beta*P.varphi;
%---------------------------------------------
j = j+1;% [Production Function] 
%---------------------------------------------
G0(j,V.c)   = -1;
G0(j,V.n)   = 1;
%---------------------------------------------
j = j+1;% [Monetary Policy] 
%---------------------------------------------
G0(j,V.r)     = 1;
G0(j,V.inf)   = -P.phi_pi;
%---------------------------------------------
j = j+1;% [Demand Shock] 
%---------------------------------------------
G0(j,V.del)   = 1;
G1(j,V.del)   = P.rho;
Psi(j,V.epsd) = P.sigma;
%---------------------------------------------
%   Solve Linear Model
%---------------------------------------------
[T,~,M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);     
