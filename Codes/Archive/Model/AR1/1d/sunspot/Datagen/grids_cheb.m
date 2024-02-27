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
G.sun_grid = linspace(O.sunbound(1),O.sunbound(2),O.sun_pts);

%--------------------------------------------------------------------------
% Weights for numerical integration from truncated normal
%--------------------------------------------------------------------------
[e_nodes,G.e_weight] = ghquad(O.e_pts);
G.e_nodes = (2^.5) * P.sigma * e_nodes;   
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
% [G.sun, G.del] = ndgrid(G.sun_grid, G.del_grid);
% G.del_nodes = numel(G.del_grid);
% G.sun_nodes = numel(G.sun_grid);
% G.nodes = numel(ndgrid(G.sun_grid, G.del_grid));
% G.griddim = size(ndgrid(G.sun_grid, G.del_grid));

G.sun = ndgrid(G.sun_grid);
G.sun_nodes = numel(G.sun_grid);
G.del = ndgrid(G.del_grid);
G.del_nodes = numel(G.del_grid);
G.nodes = numel(ndgrid(G.del_grid));
G.griddim = size(ndgrid(G.del_grid));