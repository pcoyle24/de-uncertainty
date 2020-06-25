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
pf.c_s = ones(G.griddim)*S.c_s;
pf.inf_s = ones(G.griddim)*S.inf_s;
pf.r_s = ones(G.griddim)*S.r_s;
pf.n_s = ones(G.griddim)*S.n_s;
pf.y_s = ones(G.griddim)*S.y_s;
pf.v_s = ones(G.griddim)*S.v_s;
pf.v_zlb_s = ones(G.griddim)*S.v_s;

range = 1*(chebspace(O.delbound(2),O.delbound(1),O.del_pts)'-1);

pf.c_d = ones(G.griddim)*S.c_d + range;
pf.inf_d = ones(G.griddim)*S.inf_d + range;
pf.r_d = ones(G.griddim)*S.r_d + range;
pf.n_d = ones(G.griddim)*S.n_d + range;
pf.y_d = ones(G.griddim)*S.y_d + range;
pf.v_d = ones(G.griddim)*S.v_d + range;

pf.c_zlb_d = pf.c_d;
pf.inf_zlb_d = pf.inf_d;
pf.r_zlb_d = pf.r_d;
pf.n_zlb_d = pf.n_d;
pf.y_zlb_d = pf.y_d;
pf.v_zlb_d = pf.v_d;
end


