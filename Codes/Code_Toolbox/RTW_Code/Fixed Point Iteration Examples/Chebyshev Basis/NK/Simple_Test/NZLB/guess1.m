function pf = guess1(O,P,S,G)

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
    G.del = ndgrid(G.del_grid);
    del_gr_per = (G.del - S.del) ./ S.del;
    
    % Calculate linear policy functions on discretized state space
    log_lin_pf_c   = zeros(G.griddim);
    log_lin_pf_r   = zeros(G.griddim);
    log_lin_pf_inf = zeros(G.griddim);
    
    state = [del_gr_per(:)]';
    
    log_lin_pf_c(:)   = T(V.c,[V.del])*state;
    log_lin_pf_r(:)   = T(V.r,[V.del])*state;    
    log_lin_pf_inf(:) = T(V.inf,[V.del])*state;
    
    % Convert back to levels
    pf.c    = S.c * (1 + log_lin_pf_c);
    pf.r    = S.r * (1 + log_lin_pf_r);
    pf.inf  = S.inf * (1 + log_lin_pf_inf);
    
elseif strcmp(O.loadpf,'current')
    load('pf.mat');
end