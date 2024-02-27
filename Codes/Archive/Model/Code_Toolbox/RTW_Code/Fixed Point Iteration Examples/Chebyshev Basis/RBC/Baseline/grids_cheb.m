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
% Define grid points for capital
G.k_grid = chebspace(O.kbound(1),O.kbound(2),O.k_pts);

% Define grid points for the TFP
G.z_grid = chebspace(O.zbound(1),O.zbound(2),O.z_pts);
%--------------------------------------------------------------------------
% Weights for numerical integration from truncated normal
%--------------------------------------------------------------------------
[e_nodes,G.e_weight] = ghquad(O.e_pts);
G.e_nodes = (2^.5) * sqrt(P.sig_eps) * e_nodes;    
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
[G.k_gr G.z_gr] = ndgrid(G.k_grid, G.z_grid);
G.nodes = numel(G.k_gr);
G.griddim = size(G.k_gr);