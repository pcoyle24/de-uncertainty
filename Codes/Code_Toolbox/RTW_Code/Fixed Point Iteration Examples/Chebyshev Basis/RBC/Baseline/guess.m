function pf = guess(O,P,S,G)

% pf = guess(O,P,S,G) 
%   Sets the initial policy functions
% Inputs:
%   O : structure of options
%   P : structure of parameters
%   S : structure of steady state values 
%   G : structure of grids
% Outputs:
%  pf : structure of policy functions

if strcmp(O.loadpf,'guess')
    % Retrieve log-linear solution
    V = variables;
    [T,~,eu] = linmodel(P,S,V);
    disp(['eu:    ' mat2str(eu)])
    
    % Transform discretized state space to percent deviation from steady state   
    k_gr_per = (G.k_gr - S.k) ./ S.k;  
    z_gr_per = G.z_gr - 1;

    % Calculate log-linear policy functions on discretized state space
    linpf_n = zeros(G.griddim);
    linpf_n(:) = T(V.n,[V.k V.z])*[k_gr_per(:) z_gr_per(:)]';    

    % Convert back to levels
    pf.n = P.n * (1 + linpf_n);
elseif strcmp(O.loadpf,'current')
    load 'pf.mat';
end 