function pf = guess_flat(O,P,S,G)

% pf = guess(O,P,S,G) 
%   Sets the initial policy functions (Flat line function at SS)
% Inputs:
%   O : structure of options
%   P : structure of parameters
%   S : structure of steady state values 
%   G : structure of grids
% Outputs:
%  pf : structure of policy functions
pf.c = ones(G.griddim)*S.c;
pf.inf = ones(G.griddim)*S.inf;
pf.r = ones(G.griddim)*S.r;


% pf.c = ones(G.griddim)*0.9;
% pf.inf = ones(G.griddim)*0.9;
% pf.r = ones(G.griddim)*0.9;
end