function G = grids_cheb(O,P)

% G = grids_cheb(O,P)
%   Constructs discretized state space based on the zeros of the Chebyshev
%   polynomials
% Inputs:
%     O     :   Structure of user-specified options
%     P     :   Structure of parameters
% Output:
%     G     :   Structure of grid points

%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for prefernce
G.del_grid = chebspace(O.delbound(1),O.delbound(2),O.del_pts);

%--------------------------------------------------------------------------
% Weights for numerical integration from truncated normal
%--------------------------------------------------------------------------
[e_nodes,G.e_weight] = ghquad(O.e_pts);
G.e_nodes = (2^.5) * P.sigma * e_nodes;    
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
G.del_gr = ndgrid(G.del_grid);
G.nodes = numel(G.del_gr);
G.griddim = size(G.del_gr);