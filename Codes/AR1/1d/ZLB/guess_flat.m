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

range = 0*(chebspace(O.delbound(2),O.delbound(1),O.del_pts)'-1);

pf.c = ones(G.griddim)*S.c + range;
pf.inf = ones(G.griddim)*S.inf + range;
pf.r = ones(G.griddim)*S.r + range;
pf.n = ones(G.griddim)*S.n + range;
pf.y = ones(G.griddim)*S.y + range;

% For Value Function Iteration
pf.v = ones(G.griddim)*S.v;
pf.v_zlb = ones(G.griddim)*S.v;

end