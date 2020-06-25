function S = steadystate(P)

% S = steadystate(P)
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% FOC bonds and MP Rule
    % This ultimtely may have to come above (in the code) before the section
    % Production Function and FOC Labor. It may cause errors...let's see.  
S.inf = 1;
S.r = (S.inf/P.beta);

%Production Function and FOC Labor
S.n = ((P.theta-1)/(P.theta*((1-(P.varphi/2)*(S.inf-1)^2)^P.chic)))^(1/(P.chin+P.chic)); 

S.c = S.n*(1 - (P.varphi/2)*(S.inf-1)^2);

%Demand Shock
S.del = 1; %In shock process, set \delta_t=\delta_{t-1}=\delta_{ss}; solve for \delta_{ss} 

