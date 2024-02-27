function G = grids(O,P,S)

% G = grids(O,P)
%   Constructs an evenly-spaced discretized state space
% Inputs:
%	O : Structure of user-specified options
%   P : Structure of parameters
%   S : Structure of steady state
% Output:
%   G : Structure of grid points

%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for government debt
G.del_grid = linspace(O.delbound(1),O.delbound(2),O.del_pts);

%--------------------------------------------------------------------------
% Weights for numerical integration from truncated normal
%(Don't Change) ghquad is a prewritten function
%--------------------------------------------------------------------------
[e_nodes,G.e_weight] = ghquad(O.e_pts);
G.e_nodes = (2^.5) * P.sigma * e_nodes;
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
G.del = ndgrid(G.del_grid);
G.nodes = numel(G.del);
G.griddim = size(G.del);

% Linear model -- used for guess
%G.l.del_gr = ndgrid(G.del_grid);
%G.l.nodes = numel(G.l.del_gr);
%G.l.griddim = size(G.l.del_gr);
% Non-linear model
%G.nl.del_gr = ndgrid(G.del_grid);
%G.nl.nodes = numel(G.nl.del_gr);
%G.nl.griddim = size(G.nl.del_gr);