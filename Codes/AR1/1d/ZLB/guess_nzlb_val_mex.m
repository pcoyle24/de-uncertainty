function pf = guess_nzlb_val_mex(P,S,G,O,C)
% pf = guess_nzlb(O,S,P,G,C)
% 
% Description: Calculate the policy functions of the economy absent a Zero
% Lower Bound. These will be used as the initial guesses for the economy
% with the ZLB
% Inputs:
%   P       :   Structure of parameters
%   S       :   Structure of steady state values
%   G       :   Structure of grids
%   O       :   Structure of options
%   C       :   Structure of Chebyshev polynomial parameters
% Output:
%   pf      :   Structure of Policy Functions  


global pi_yesterday

% Retrieve initial policy functions
pf = guess_flat(O,P,S,G); 

%--------------------------------------------------------------------------
% Initialize algorithm
%--------------------------------------------------------------------------
% Calculate initial Chebyshev coefficients          
C.Ac = Fchebweights11(O.n1,O.del_pts,pf.c,C.T,C.P,C.X);
C.Ainf = Fchebweights11(O.n1,O.del_pts,pf.inf,C.T,C.P,C.X);
C.Av = Fchebweights11(O.n1,O.del_pts,pf.v,C.T,C.P,C.X);

theta_old = [C.Ac,C.Ainf,C.Av];


options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',true);
display(char('Calculating the Policy Functions for Economy without ZLB'))
func    = @(weights) eqm_nzlb_val_mex(weights,P,S,G,O,C);

[theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func,theta_old,[],[],options);

C.Ac = theta_new(:,1);
C.Ainf = theta_new(:,2);
C.Av = theta_new(:,3);

pf.c        = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ac,C.max,C.T,C.P);     
pf.inf      = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ainf,C.max,C.T,C.P); 
pf.pitilde  = pf.inf/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
pf.n        = (pf.c./(1-(P.varphi/2).*(pf.pitilde-1).^2));
pf.y        = pf.n;
pf.w        = pf.n.^P.chin.*pf.c.^P.chic;
pf.r        = S.pi_targ/P.beta*((pf.inf./S.pi_targ).^(P.phi_pi).*(pf.y./S.y).^(P.phi_y));
pf.v        = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av,C.max,C.T,C.P); 

pf.c_zlb        = pf.c;
pf.inf_zlb      = pf.inf;
pf.pitilde_zlb  = pf.pitilde;
pf.n_zlb        = pf.n;
pf.y_zlb        = pf.y;
pf.w_zlb        = pf.w;
pf.r_zlb        = ones(size(pf.c));
pf.v_zlb        = pf.v;

display(char('Policy Functions Calculated'))
    

    
    
