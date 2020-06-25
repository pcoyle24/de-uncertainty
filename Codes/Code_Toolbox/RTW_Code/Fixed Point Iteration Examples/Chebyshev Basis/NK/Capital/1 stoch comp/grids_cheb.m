function G = grids_cheb(O,P,S)

% G = grids_cheb(O,P)
%   Constructs discretized state space based on the zeros of the Chebyshev
%   polynomials
% Inputs:
%     O     :   Structure of user-specified options
%     P     :   Structure of parameters
%     S     :   Structure of steady state values
% Output:
%     G     :   Structure of grid points

%--------------------------------------------------------------------------
% Set up grids by discretizing the state space
%--------------------------------------------------------------------------
% Define grid points for government liabilities
G.a_grid = chebspace(O.abound(1),O.abound(2),O.a_pts);

% Define grid points for government debt
G.b_grid = chebspace(O.bbound(1),O.bbound(2),O.b_pts);

% Define grid points for capital
G.k_grid = chebspace(O.kbound(1),O.kbound(2),O.k_pts);

% Define grid points for shock in tax rule
G.fp_grid = chebspace(O.fpbound(1),O.fpbound(2),O.fp_pts);
%--------------------------------------------------------------------------
% Weights for numerical integration from truncated normal
%--------------------------------------------------------------------------
[e_nodes,G.e_weight] = ghquad(O.e_pts);
G.e_nodes = (2^.5) * P.sig_fp * e_nodes;
%--------------------------------------------------------------------------
% Set up remaining grids for linear model (needed for initial guess)
%--------------------------------------------------------------------------
% Define grid points for r
rmin = S.r - 0.01;
rmax = S.r + 0.01;
G.r_pts = O.a_pts;
G.r_grid = linspace(rmin,rmax,O.a_pts);

%Set m so that it corresponds to the way a is discretized
G.m_grid = G.a_grid - G.r_grid.*G.b_grid;
%--------------------------------------------------------------------------
% Construct grid arrays
%--------------------------------------------------------------------------
% Linear model -- used for guess
[G.l.b_gr G.l.r_gr G.l.m_gr G.l.k_gr] = ...
        ndgrid(G.b_grid, G.r_grid, G.m_grid, G.k_grid);
G.l.nodes = numel(G.l.b_gr);
G.l.griddim = size(G.l.b_gr);
% Non-linear model
[G.nl.a_gr G.nl.b_gr G.nl.k_gr G.nl.fp_gr] = ...
        ndgrid(G.a_grid, G.b_grid, G.k_grid, G.fp_grid);
G.nl.nodes = numel(G.nl.b_gr);
G.nl.griddim = size(G.nl.b_gr);