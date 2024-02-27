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
C.Ac_s = Fchebweights11(O.n1,O.del_pts,pf.c_s,C.T,C.P,C.X);
C.Ainf_s = Fchebweights11(O.n1,O.del_pts,pf.inf_s,C.T,C.P,C.X);
C.Av_s = Fchebweights11(O.n1,O.del_pts,pf.v_s,C.T,C.P,C.X);

theta_old_s = [C.Ac_s,C.Ainf_s,C.Av_s];


options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',true);
display(char('Calculating the Policy Functions for Economy without ZLB'))
func    = @(weights) eqm_nzlb_val_mex(weights,P,S,G,O,C);

[theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func,theta_old_s,[],[],options);

C.Ac_s = theta_new(:,1);
C.Ainf_s = theta_new(:,2);
C.Av_s = theta_new(:,3);

pf.c_s        = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ac_s,C.max,C.T,C.P);     
pf.inf_s      = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Ainf_s,C.max,C.T,C.P); 
pf.pitilde_s  = pf.inf_s/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
pf.n_s        = (pf.c_s./(1-(P.varphi/2).*(pf.pitilde_s-1).^2));
pf.y_s        = pf.n_s;
pf.w_s        = pf.n_s.^P.chin.*pf.c_s.^P.chic;
pf.r_s        = S.pi_targ/P.beta*((pf.inf_s./S.pi_targ).^(P.phi_pi).*(pf.y_s./S.y).^(P.phi_y));
pf.v_s        = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av_s,C.max,C.T,C.P); 

pf.c_zlb_s        = pf.c_s;
pf.inf_zlb_s      = pf.inf_s;
pf.pitilde_zlb_s  = pf.pitilde_s;
pf.n_zlb_s        = pf.n_s;
pf.y_zlb_s        = pf.y_s;
pf.w_zlb_s        = pf.w_s;
pf.r_zlb_s        = ones(size(pf.c_s));
pf.v_zlb_s        = pf.v_s;

display(char('Policy Functions Calculated'))
    

    
    
