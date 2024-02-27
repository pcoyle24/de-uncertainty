function [R,J] = eqm_mex(weights,P,S,G,O,C,CondProb)

% [R,J] = eqm_jacob(weights,P,S,G,O,C,CondProb)
%   Outputs residuals and Analytical Jacobian of the equilibrium system of 
%   equations for least squares minimization/Chebyshev interpolation method.
% Inputs:
%   weights :   Parameter Weights for Cheb Poly (Matrix)
%   P       :   Structure of parameters
%   S       :   Structure of steady state values
%   G       :   Structure of grids
%   O       :   Structure of options
%   C       :   Structure of Chebyshev polynomial parameters
%   CondProb:   Sun Spot Transition Matrix
% Output:
%   R       :   Residuals
%   J       :   Jacobian

global pi_yesterday

% Allocate Space for Residuals

res1_s        = zeros(G.nodes,1);
res2_s        = zeros(G.nodes,1);
res1_zlb_s        = zeros(G.nodes,1);
res2_zlb_s        = zeros(G.nodes,1);
res1_d        = zeros(G.nodes,1);
res2_d        = zeros(G.nodes,1);
res1_zlb_d        = zeros(G.nodes,1);
res2_zlb_d        = zeros(G.nodes,1);

[~,col] = size(weights);

lhs_jacob_ee_s = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_ee_d = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_ee_s = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_ee_d = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_pc_s = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_pc_d = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_pc_s = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_pc_d = zeros(G.nodes,col*(O.n1+1));


lhs_jacob_ee_zlb_s = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_ee_zlb_d = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_ee_zlb_s = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_ee_zlb_d = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_pc_zlb_s = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_pc_zlb_d = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_pc_zlb_s = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_pc_zlb_d = zeros(G.nodes,col*(O.n1+1));


% Map Coefficients 
C_Ac_s    = weights(1:O.n1+1,1);
C_Ainf_s  = weights(1:O.n1+1,2);
C_Ac_zlb_s    = weights(1:O.n1+1,3);
C_Ainf_zlb_s  = weights(1:O.n1+1,4);

C_Ac_d    = weights(1:O.n1+1,5);
C_Ainf_d  = weights(1:O.n1+1,6);
C_Ac_zlb_d    = weights(1:O.n1+1,7);
C_Ainf_zlb_d  = weights(1:O.n1+1,8);

% Get Partial Derivatives for Analytical Jacobian
get_partial_derivs;

% Get Policy Functions
% Standard Regime, NZLB
 pf_c_s           = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ac_s,C.max,C.T,C.P);    
 pf_inf_s         = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ainf_s,C.max,C.T,C.P);
 pf_pitilde_s  = pf_inf_s/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
 pf_n_s        = (pf_c_s./(1-(P.varphi/2).*(pf_pitilde_s-1).^2));
 pf_y_s        = pf_n_s;
 pf_w_s        = pf_n_s.^P.chin.*pf_c_s.^P.chic;
 pf_r_s        = S.pi_targ/P.beta*((pf_inf_s./S.pi_targ).^(P.phi_pi).*(pf_y_s./S.y).^(P.phi_y));
 
% Standard Regime, ZLB
 pf_c_zlb_s       = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ac_zlb_s,C.max,C.T,C.P);    
 pf_inf_zlb_s     = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ainf_zlb_s,C.max,C.T,C.P);
 pf_r_zlb_s       = ones(size(pf_r_s));
 
% Deflationary Regime, NZLB
 pf_c_d           = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ac_d,C.max,C.T,C.P);    
 pf_inf_d         = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ainf_d,C.max,C.T,C.P);
 pf_pitilde_d  = pf_inf_d/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
 pf_n_d        = (pf_c_d./(1-(P.varphi/2).*(pf_pitilde_d-1).^2));
 pf_y_d        = pf_n_d;
 pf_w_d        = pf_n_d.^P.chin.*pf_c_d.^P.chic;
 pf_r_d        = S.pi_targ/P.beta*((pf_inf_d./S.pi_targ).^(P.phi_pi).*(pf_y_d./S.y).^(P.phi_y));

 % Deflationary Regime, ZLB
 pf_c_zlb_d       = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ac_zlb_d,C.max,C.T,C.P);    
 pf_inf_zlb_d     = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C_Ainf_zlb_d,C.max,C.T,C.P);
 pf_r_zlb_d       = ones(size(pf_r_s));
 
 % Get r_today's cheb weights
 C_Ar_s = Fchebweights11(O.n1,O.del_pts,pf_r_s,C.T,C.P,C.X);
 C_Ar_d = Fchebweights11(O.n1,O.del_pts,pf_r_d,C.T,C.P,C.X);

 for i = 1:G.nodes
     del_today         = G.del(i);
     for k = 1:G.sun_nodes
         cond_prob_sun = CondProb(k,:);
         if k == 1 %Standard Regime
             c_temp            = pf_c_s(i);
             inf_temp          = pf_inf_s(i);

             pi_tilde_temp     = inf_temp/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
             n_temp            = (c_temp/(1-(P.varphi/2)*(pi_tilde_temp-1)^2));
             y_temp            = n_temp;
             w_temp            = n_temp^P.chin*c_temp^P.chic;
             r_temp            = S.pi_targ/P.beta*((inf_temp/S.pi_targ)^(P.phi_pi)*(y_temp/S.y)^(P.phi_y));
             if r_temp  >= 1
                 c_today        = c_temp;
                 inf_today      = inf_temp;

                 pi_tilde_today = pi_tilde_temp;
                 n_today        = n_temp;
                 y_today        = y_temp;
                 w_today        = w_temp;
                 r_today        = r_temp;

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                 %For GH Integration
                exp_ee_s_int_temp = 0;
                exp_pc_s_int_temp = 0;
                exp_ee_d_int_temp = 0;
                exp_pc_d_int_temp = 0;
                exp_ee_s_int = 0;
                exp_pc_s_int = 0;

                dee_s_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_s_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_s_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_s_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_s_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_s_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_s_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_s_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                dpc_s_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_s_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_s_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_s_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_s_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_s_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_s_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_s_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_s,C.max,C.T,C.P);
                                 inf_s_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow     = inf_s_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow            = (c_s_tomorrow/(1-(P.varphi/2)*(pi_tilde_s_tomorrow-1)^2));
                                 y_s_tomorrow            = n_s_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_s             = c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-1);
                                 exp_pc_s             = (y_s_tomorrow/c_s_tomorrow^(P.chic))*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow;
                                 exp_ee_s_jacob = exp_ee_s;


                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_s_dax_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow*...
                                                (c_s_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_s_tomorrow) - P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_s_tomorrow);

                                 % Derivatives with respect to b coefficients
                                 dpc_s_dbx_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_s_tomorrow,pi_tilde_s_tomorrow)*(pi_tilde_s_tomorrow-1)*pi_tilde_s_tomorrow) + y_s_tomorrow*(2*pi_tilde_s_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));

                                 exp_ee_s_int_temp         = exp_ee_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_s);
                                 exp_pc_s_int_temp         = exp_pc_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_s);

                                 dpc_s_dax_s_int = dpc_s_dax_s_int + pi^(-0.5)*G.e_weight(j)*dpc_s_dax_s;
                                 dpc_s_dbx_s_int = dpc_s_dbx_s_int + pi^(-0.5)*G.e_weight(j)*dpc_s_dbx_s;                 

                             else 
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_s,C.max,C.T,C.P);
                                 inf_s_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow_zlb     = inf_s_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow_zlb            = (c_s_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_s_tomorrow_zlb-1)^2));
                                 y_s_tomorrow_zlb            = n_s_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_s             = c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-1);
                                 exp_pc_zlb_s             = (y_s_tomorrow_zlb/c_s_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb;
                                 exp_ee_s_jacob = exp_ee_zlb_s;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dee_s_daxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today*P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_s_tomorrow_zlb^(-1);

                                 % Derivatives with respect to b coefficients
                                 dee_s_dbxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today*c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_s_daxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb*...
                                                (c_s_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_s_tomorrow_zlb) - P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_s_tomorrow_zlb);

                                 % Derivatives with respect to b coefficients
                                 dpc_s_dbxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_s_tomorrow_zlb,pi_tilde_s_tomorrow_zlb)*(pi_tilde_s_tomorrow_zlb-1)*pi_tilde_s_tomorrow_zlb) + y_s_tomorrow_zlb*(2*pi_tilde_s_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                    

                                 exp_ee_s_int_temp         = exp_ee_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_s);
                                 exp_pc_s_int_temp         = exp_pc_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_s);

                                 dee_s_daxzlb_s_int = dee_s_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dee_s_daxzlb_s;
                                 dee_s_dbxzlb_s_int = dee_s_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dee_s_dbxzlb_s;

                                 dpc_s_daxzlb_s_int = dpc_s_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpc_s_daxzlb_s;
                                 dpc_s_dbxzlb_s_int = dpc_s_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpc_s_dbxzlb_s;
                             end 
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_d,C.T,C.P);

                                 % Get Interpolated Values
                                 %Deflationary Regime
                                 c_d_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_d,C.max,C.T,C.P);
                                 inf_d_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow     = inf_d_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow            = (c_d_tomorrow/(1-(P.varphi/2)*(pi_tilde_d_tomorrow-1)^2));
                                 y_d_tomorrow            = n_d_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_d             = c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-1);
                                 exp_pc_d             = (y_d_tomorrow/c_d_tomorrow^(P.chic))*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow;
                                 exp_ee_d_jacob = exp_ee_d;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dee_s_dax_d = -cond_prob_sun(2)*P.beta*del_today*r_today*P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_d_tomorrow^(-1);                     

                                 % Derivatives with respect to b coefficients
                                 dee_s_dbx_d = -cond_prob_sun(2)*P.beta*del_today*r_today*c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-2)*dpi_dbx(basis_interp);   

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_s_dax_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow*...
                                                (c_d_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_d_tomorrow) - P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_d_tomorrow);           

                                 % Derivatives with respect to b coefficients
                                 dpc_s_dbx_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_d_tomorrow,pi_tilde_d_tomorrow)*(pi_tilde_d_tomorrow-1)*pi_tilde_d_tomorrow) + y_d_tomorrow*(2*pi_tilde_d_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));

                                 exp_ee_d_int_temp         = exp_ee_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_d);
                                 exp_pc_d_int_temp         = exp_pc_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_d);

                                 dee_s_dax_d_int = dee_s_dax_d_int + pi^(-0.5)*G.e_weight(j)*dee_s_dax_d;
                                 dee_s_dbx_d_int = dee_s_dbx_d_int + pi^(-0.5)*G.e_weight(j)*dee_s_dbx_d;

                                 dpc_s_dax_d_int = dpc_s_dax_d_int + pi^(-0.5)*G.e_weight(j)*dpc_s_dax_d;
                                 dpc_s_dbx_d_int = dpc_s_dbx_d_int + pi^(-0.5)*G.e_weight(j)*dpc_s_dbx_d;
                             else
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_d,C.T,C.P);

                                 %Deflationary Regime
                                 c_d_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_d,C.max,C.T,C.P);
                                 inf_d_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow_zlb     = inf_d_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow_zlb            = (c_d_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_d_tomorrow_zlb-1)^2));
                                 y_d_tomorrow_zlb            = n_d_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_d             = c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-1);
                                 exp_pc_zlb_d             = (y_d_tomorrow_zlb/c_d_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb;
                                 exp_ee_d_jacob = exp_ee_zlb_d;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dee_s_daxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today*P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_d_tomorrow_zlb^(-1);                     

                                 % Derivatives with respect to b coefficients
                                 dee_s_dbxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today*c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);   

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_s_daxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb*...
                                                (c_d_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_d_tomorrow_zlb) - P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_d_tomorrow_zlb);           

                                 % Derivatives with respect to b coefficients
                                 dpc_s_dbxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_d_tomorrow_zlb,pi_tilde_d_tomorrow_zlb)*(pi_tilde_d_tomorrow_zlb-1)*pi_tilde_d_tomorrow_zlb) + y_d_tomorrow_zlb*(2*pi_tilde_d_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                               


                                 exp_ee_d_int_temp         = exp_ee_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_d);
                                 exp_pc_d_int_temp         = exp_pc_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_d);

                                 dee_s_daxzlb_d_int = dee_s_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dee_s_daxzlb_d;
                                 dee_s_dbxzlb_d_int = dee_s_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dee_s_dbxzlb_d;

                                 dpc_s_daxzlb_d_int = dpc_s_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpc_s_daxzlb_d;
                                 dpc_s_dbxzlb_d_int = dpc_s_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpc_s_dbxzlb_d; 

                             end              
                         end
                     end
                     exp_ee_s_int = cond_prob_sun*[exp_ee_s_int_temp;exp_ee_d_int_temp];
                     exp_pc_s_int = cond_prob_sun*[exp_pc_s_int_temp;exp_pc_d_int_temp];

                     if r_s_tomorrow >= 1
                         % RHS derivative of EE (for Jacobian)
                         % Derivatives with respect to a coefficients
                         dee_s_dax_s = P.beta*del_today*(dr_dax(C.basis(:,i),inf_today,y_today,pi_tilde_today)*(cond_prob_sun*[exp_ee_s_jacob;exp_ee_d_jacob])...
                                     - cond_prob_sun(1)*r_today*P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_s_tomorrow^(-1));

                         % Derivatives with respect to b coefficients
                         dee_s_dbx_s = P.beta*del_today*(dr_dbx(C.basis(:,i),inf_today,y_today,c_today,pi_tilde_today)*(cond_prob_sun*[exp_ee_s_jacob;exp_ee_d_jacob])...
                                     - cond_prob_sun(1)*r_today*c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-2)*dpi_dbx(basis_interp));

                         dee_s_dax_s_int = dee_s_dax_s_int + pi^(-0.5)*G.e_weight(j)*dee_s_dax_s;
                         dee_s_dbx_s_int = dee_s_dbx_s_int + pi^(-0.5)*G.e_weight(j)*dee_s_dbx_s;
                     end
                 end
                 % Residuals
                 res1_s(i) = c_today^(-P.chic) - P.beta*del_today*r_today*exp_ee_s_int;
                 res2_s(i) = (y_today/c_today^(P.chic))*(P.varphi*(pi_tilde_today-1)*pi_tilde_today-(1-P.theta)-P.theta*(1-P.tau)*w_today) - P.beta*del_today*exp_pc_s_int;
                 
                 % Jacobian
                 % Jacobian enteries relating to EE
                 lhs_jacob_ee_s(i,:) = [-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i)), 0*-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))...
                                    0*-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i)), 0*-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))];
                 rhs_jacob_ee_s(i,:) = [dee_s_dax_s_int, dee_s_dbx_s_int, dee_s_daxzlb_s_int, dee_s_dbxzlb_s_int ...
                                      dee_s_dax_d_int, dee_s_dbx_d_int, dee_s_daxzlb_d_int, dee_s_dbxzlb_d_int];

                 % Jacobian enteries relating to PC
                 lhs_jacob_pc_ax = (dy_dax(C.basis(:,i),pi_tilde_today)*c_today^(-P.chic) - P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i))*y_today)* ...
                                  (P.varphi*(pi_tilde_today - 1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
                 lhs_jacob_pc_ax = lhs_jacob_pc_ax + -P.theta*(1-P.tau)*dw_dax(C.basis(:,i),c_today,pi_tilde_today,n_today)* ...
                                  y_today/(c_today^(P.chic));

                 lhs_jacob_pc_bx = dy_dbx(C.basis(:,i),c_today,pi_tilde_today)/c_today^(P.chic)*(P.varphi*(pi_tilde_today-1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
                 lhs_jacob_pc_bx = lhs_jacob_pc_bx + y_today/(c_today^(P.chic))*(P.varphi*(2*pi_tilde_today*dpitilde_dbx(C.basis(:,i))-dpitilde_dbx(C.basis(:,i))) - P.theta*(1-P.tau)*dw_dbx(C.basis(:,i),c_today,pi_tilde_today,n_today));

                 lhs_jacob_pc_s(i,:)  = [lhs_jacob_pc_ax, lhs_jacob_pc_bx, 0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx...
                                       0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx, 0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx];
                 rhs_jacob_pc_s(i,:)  = [dpc_s_dax_s_int,dpc_s_dbx_s_int, dpc_s_daxzlb_s_int,dpc_s_dbxzlb_s_int...
                                         dpc_s_dax_d_int,dpc_s_dbx_d_int, dpc_s_daxzlb_d_int,dpc_s_dbxzlb_d_int];                
             else
                 c_today_zlb        = pf_c_zlb_s(i);
                 inf_today_zlb      = pf_inf_zlb_s(i);

                 pi_tilde_today_zlb     = inf_today_zlb/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
                 n_today_zlb           = (c_today_zlb/(1-(P.varphi/2)*(pi_tilde_today_zlb-1)^2));
                 y_today_zlb            = n_today_zlb;
                 w_today_zlb            = n_today_zlb^P.chin*c_today_zlb^P.chic;
                 r_today_zlb            = S.r_zlb;

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                %For GH Integration
                exp_ee_s_int_zlb_temp = 0;
                exp_pc_s_int_zlb_temp = 0;
                exp_ee_d_int_zlb_temp = 0;
                exp_pc_d_int_zlb_temp = 0;
                exp_ee_s_int_zlb = 0;
                exp_pc_s_int_zlb = 0;

                deezlb_s_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_s_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_s_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_s_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_s_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_s_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_s_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_s_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                dpczlb_s_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_s_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_s_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_s_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_s_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_s_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_s_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_s_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_s,C.max,C.T,C.P);
                                 inf_s_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow     = inf_s_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow            = (c_s_tomorrow/(1-(P.varphi/2)*(pi_tilde_s_tomorrow-1)^2));
                                 y_s_tomorrow            = n_s_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_s             = c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-1);
                                 exp_pc_s             = (y_s_tomorrow/c_s_tomorrow^(P.chic))*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_s_dax_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_s_tomorrow^(-1);

                                 % Derivatives with respect to b coefficients
                                 deezlb_s_dbx_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-2)*dpi_dbx(basis_interp);

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_s_dax_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow*...
                                                (c_s_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_s_tomorrow) - P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_s_tomorrow);

                                % Derivatives with respect to b coefficients
                                 dpczlb_s_dbx_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_s_tomorrow,pi_tilde_s_tomorrow)*(pi_tilde_s_tomorrow-1)*pi_tilde_s_tomorrow) + y_s_tomorrow*(2*pi_tilde_s_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));                    

                                 exp_ee_s_int_zlb_temp         = exp_ee_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_s);
                                 exp_pc_s_int_zlb_temp         = exp_pc_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_s);

                                 deezlb_s_dax_s_int = deezlb_s_dax_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_dax_s;
                                 deezlb_s_dbx_s_int = deezlb_s_dbx_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_dbx_s;

                                 dpczlb_s_dax_s_int = dpczlb_s_dax_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_dax_s;
                                 dpczlb_s_dbx_s_int = dpczlb_s_dbx_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_dbx_s;

                             else 
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_s,C.max,C.T,C.P);
                                 inf_s_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow_zlb     = inf_s_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow_zlb            = (c_s_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_s_tomorrow_zlb-1)^2));
                                 y_s_tomorrow_zlb            = n_s_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_s             = c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-1);
                                 exp_pc_zlb_s             = (y_s_tomorrow_zlb/c_s_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_s_daxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_s_tomorrow_zlb^(-1);

                                 % Derivatives with respect to b coefficients
                                 deezlb_s_dbxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_s_daxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb*...
                                                (c_s_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_s_tomorrow_zlb) - P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_s_tomorrow_zlb);

                                 % Derivatives with respect to b coefficients
                                 dpczlb_s_dbxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_s_tomorrow_zlb,pi_tilde_s_tomorrow_zlb)*(pi_tilde_s_tomorrow_zlb-1)*pi_tilde_s_tomorrow_zlb) + y_s_tomorrow_zlb*(2*pi_tilde_s_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                    


                                 exp_ee_s_int_zlb_temp         = exp_ee_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_s);
                                 exp_pc_s_int_zlb_temp         = exp_pc_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_s);

                                 deezlb_s_daxzlb_s_int = deezlb_s_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_daxzlb_s;
                                 deezlb_s_dbxzlb_s_int = deezlb_s_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_dbxzlb_s;

                                 dpczlb_s_daxzlb_s_int = dpczlb_s_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_daxzlb_s;
                                 dpczlb_s_dbxzlb_s_int = dpczlb_s_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_dbxzlb_s;
                             end 
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_d,C.T,C.P);

                                 % Get Interpolated Values
                                 %Deflationary Regime
                                 c_d_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_d,C.max,C.T,C.P);
                                 inf_d_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow     = inf_d_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow            = (c_d_tomorrow/(1-(P.varphi/2)*(pi_tilde_d_tomorrow-1)^2));
                                 y_d_tomorrow            = n_d_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_d             = c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-1);
                                 exp_pc_d             = (y_d_tomorrow/c_d_tomorrow^(P.chic))*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_s_dax_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_d_tomorrow^(-1);                     

                                 % Derivatives with respect to b coefficients
                                 deezlb_s_dbx_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-2)*dpi_dbx(basis_interp);   

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_s_dax_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow*...
                                                (c_d_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_d_tomorrow) - P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_d_tomorrow);           

                                 % Derivatives with respect to b coefficients
                                 dpczlb_s_dbx_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_d_tomorrow,pi_tilde_d_tomorrow)*(pi_tilde_d_tomorrow-1)*pi_tilde_d_tomorrow) + y_d_tomorrow*(2*pi_tilde_d_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));                               

                                 exp_ee_d_int_zlb_temp         = exp_ee_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_d);
                                 exp_pc_d_int_zlb_temp         = exp_pc_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_d);

                                 deezlb_s_dax_d_int = deezlb_s_dax_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_dax_d;
                                 deezlb_s_dbx_d_int = deezlb_s_dbx_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_dbx_d;

                                 dpczlb_s_dax_d_int = dpczlb_s_dax_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_dax_d;
                                 dpczlb_s_dbx_d_int = dpczlb_s_dbx_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_dbx_d;
                             else
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_d,C.T,C.P);

                                 %Deflationary Regime
                                 c_d_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_d,C.max,C.T,C.P);
                                 inf_d_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow_zlb     = inf_d_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow_zlb            = (c_d_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_d_tomorrow_zlb-1)^2));
                                 y_d_tomorrow_zlb            = n_d_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_d             = c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-1);
                                 exp_pc_zlb_d             = (y_d_tomorrow_zlb/c_d_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_s_daxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_d_tomorrow_zlb^(-1);                     

                                 % Derivatives with respect to b coefficients
                                 deezlb_s_dbxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);   

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_s_daxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb*...
                                                (c_d_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_d_tomorrow_zlb) - P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_d_tomorrow_zlb);           

                                 % Derivatives with respect to b coefficients
                                 dpczlb_s_dbxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_d_tomorrow_zlb,pi_tilde_d_tomorrow_zlb)*(pi_tilde_d_tomorrow_zlb-1)*pi_tilde_d_tomorrow_zlb) + y_d_tomorrow_zlb*(2*pi_tilde_d_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                               

                                 exp_ee_d_int_zlb_temp         = exp_ee_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_d);
                                 exp_pc_d_int_zlb_temp         = exp_pc_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_d);

                                 deezlb_s_daxzlb_d_int = deezlb_s_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_daxzlb_d;
                                 deezlb_s_dbxzlb_d_int = deezlb_s_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_s_dbxzlb_d;

                                 dpczlb_s_daxzlb_d_int = dpczlb_s_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_daxzlb_d;
                                 dpczlb_s_dbxzlb_d_int = dpczlb_s_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_s_dbxzlb_d; 

                             end              
                         end
                     end
                     exp_ee_s_int_zlb = cond_prob_sun*[exp_ee_s_int_zlb_temp;exp_ee_d_int_zlb_temp];
                     exp_pc_s_int_zlb = cond_prob_sun*[exp_pc_s_int_zlb_temp;exp_pc_d_int_zlb_temp];         
                 end

                 % Residuals
                 res1_zlb_s(i) = c_today_zlb^(-P.chic) - P.beta*del_today*r_today_zlb*exp_ee_s_int_zlb;
                 res2_zlb_s(i) = (y_today_zlb/c_today_zlb^(P.chic))*(P.varphi*(pi_tilde_today_zlb-1)*pi_tilde_today_zlb-(1-P.theta)-P.theta*(1-P.tau)*w_today_zlb) - P.beta*del_today*exp_pc_s_int_zlb;
                 
                 % Jacobian
                 % Jacobian enteries relating to EE
                 lhs_jacob_ee_zlb_s(i,:) = [0*-P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i)), -P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i))...
                                    0*-P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i)), 0*-P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))];
                 rhs_jacob_ee_zlb_s(i,:) = [deezlb_s_dax_s_int, deezlb_s_dbx_s_int, deezlb_s_daxzlb_s_int, deezlb_s_dbxzlb_s_int ...
                                      deezlb_s_dax_d_int, deezlb_s_dbx_d_int, deezlb_s_daxzlb_d_int, deezlb_s_dbxzlb_d_int];

                 % Jacobian enteries relating to PC
                 lhs_jacob_pc_ax = (dy_dax(C.basis(:,i),pi_tilde_today_zlb)*c_today_zlb^(-P.chic) - P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i))*y_today_zlb)* ...
                                  (P.varphi*(pi_tilde_today_zlb - 1)*pi_tilde_today_zlb - (1-P.theta) - P.theta*(1-P.tau)*w_today_zlb);
                 lhs_jacob_pc_ax = lhs_jacob_pc_ax + -P.theta*(1-P.tau)*dw_dax(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb,n_today_zlb)* ...
                                  y_today_zlb/(c_today_zlb^(P.chic));

                 lhs_jacob_pc_bx = dy_dbx(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb)/c_today_zlb^(P.chic)*(P.varphi*(pi_tilde_today_zlb-1)*pi_tilde_today_zlb - (1-P.theta) - P.theta*(1-P.tau)*w_today_zlb);
                 lhs_jacob_pc_bx = lhs_jacob_pc_bx + y_today_zlb/(c_today_zlb^(P.chic))*(P.varphi*(2*pi_tilde_today_zlb*dpitilde_dbx(C.basis(:,i))-dpitilde_dbx(C.basis(:,i))) - P.theta*(1-P.tau)*dw_dbx(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb,n_today_zlb));

                 lhs_jacob_pc_zlb_s(i,:)  = [0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx, lhs_jacob_pc_ax, lhs_jacob_pc_bx...
                                             0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx, 0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx];
                 rhs_jacob_pc_zlb_s(i,:)  = [dpczlb_s_dax_s_int,dpczlb_s_dbx_s_int, dpczlb_s_daxzlb_s_int,dpczlb_s_dbxzlb_s_int...
                                             dpczlb_s_dax_d_int,dpczlb_s_dbx_d_int, dpczlb_s_daxzlb_d_int,dpczlb_s_dbxzlb_d_int];
             end
         else %Deflationary Regime
             c_temp            = pf_c_d(i);
             inf_temp          = pf_inf_d(i);

             pi_tilde_temp     = inf_temp/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
             n_temp            = (c_temp/(1-(P.varphi/2)*(pi_tilde_temp-1)^2));
             y_temp            = n_temp;
             w_temp            = n_temp^P.chin*c_temp^P.chic;
             r_temp            = S.pi_targ/P.beta*((inf_temp/S.pi_targ)^(P.phi_pi)*(y_temp/S.y)^(P.phi_y));
             if r_temp  >= 1
                 c_today        = c_temp;
                 inf_today      = inf_temp;

                 pi_tilde_today = pi_tilde_temp;
                 n_today        = n_temp;
                 y_today        = y_temp;
                 w_today        = w_temp;
                 r_today        = r_temp;

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                 %For GH Integration
                exp_ee_s_int_temp = 0;
                exp_pc_s_int_temp = 0;
                exp_ee_d_int_temp = 0;
                exp_pc_d_int_temp = 0;
                exp_ee_d_int = 0;
                exp_pc_d_int = 0;

                dee_d_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_d_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_d_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_d_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_d_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_d_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_d_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dee_d_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                dpc_d_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_d_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_d_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_d_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_d_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_d_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_d_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpc_d_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_s,C.max,C.T,C.P);
                                 inf_s_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow     = inf_s_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow            = (c_s_tomorrow/(1-(P.varphi/2)*(pi_tilde_s_tomorrow-1)^2));
                                 y_s_tomorrow            = n_s_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_s             = c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-1);
                                 exp_pc_s             = (y_s_tomorrow/c_s_tomorrow^(P.chic))*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow;
                                 exp_ee_s_jacob = exp_ee_s;


                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dee_d_dax_s = -cond_prob_sun(1)*P.beta*del_today*r_today*P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_s_tomorrow^(-1);

                                 % Derivatives with respect to b coefficients
                                 dee_d_dbx_d = -cond_prob_sun(1)*P.beta*del_today*r_today*c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-2)*dpi_dbx(basis_interp);

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_d_dax_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow*...
                                                (c_s_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_s_tomorrow) - P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_s_tomorrow);

                                 % Derivatives with respect to b coefficients
                                 dpc_d_dbx_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_s_tomorrow,pi_tilde_s_tomorrow)*(pi_tilde_s_tomorrow-1)*pi_tilde_s_tomorrow) + y_s_tomorrow*(2*pi_tilde_s_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));                    


                                 exp_ee_s_int_temp         = exp_ee_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_s);
                                 exp_pc_s_int_temp         = exp_pc_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_s);

                                 dee_d_dax_s_int = dee_d_dax_s_int + pi^(-0.5)*G.e_weight(j)*dee_d_dax_s;
                                 dee_d_dbx_s_int = dee_d_dbx_s_int + pi^(-0.5)*G.e_weight(j)*dee_d_dbx_d;

                                 dpc_d_dax_s_int = dpc_d_dax_s_int + pi^(-0.5)*G.e_weight(j)*dpc_d_dax_s;
                                 dpc_d_dbx_s_int = dpc_d_dbx_s_int + pi^(-0.5)*G.e_weight(j)*dpc_d_dbx_s;                 

                             else 
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_s,C.max,C.T,C.P);
                                 inf_s_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow_zlb     = inf_s_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow_zlb            = (c_s_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_s_tomorrow_zlb-1)^2));
                                 y_s_tomorrow_zlb            = n_s_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_s             = c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-1);
                                 exp_pc_zlb_s             = (y_s_tomorrow_zlb/c_s_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb;
                                 exp_ee_s_jacob = exp_ee_zlb_s;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dee_d_daxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today*P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_s_tomorrow_zlb^(-1);

                                 % Derivatives with respect to b coefficients
                                 dee_d_dbxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today*c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_d_daxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb*...
                                                (c_s_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_s_tomorrow_zlb) - P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_s_tomorrow_zlb);

                                 % Derivatives with respect to b coefficients
                                 dpc_d_dbxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_s_tomorrow_zlb,pi_tilde_s_tomorrow_zlb)*(pi_tilde_s_tomorrow_zlb-1)*pi_tilde_s_tomorrow_zlb) + y_s_tomorrow_zlb*(2*pi_tilde_s_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                    

                                 exp_ee_s_int_temp         = exp_ee_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_s);
                                 exp_pc_s_int_temp         = exp_pc_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_s);

                                 dee_d_daxzlb_s_int = dee_d_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dee_d_daxzlb_s;
                                 dee_d_dbxzlb_s_int = dee_d_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dee_d_dbxzlb_s;

                                 dpc_d_daxzlb_s_int = dpc_d_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpc_d_daxzlb_s;
                                 dpc_d_dbxzlb_s_int = dpc_d_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpc_d_dbxzlb_s;
                             end 
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_d,C.T,C.P);

                                 % Get Interpolated Values
                                 %Deflationary Regime
                                 c_d_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_d,C.max,C.T,C.P);
                                 inf_d_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow     = inf_d_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow            = (c_d_tomorrow/(1-(P.varphi/2)*(pi_tilde_d_tomorrow-1)^2));
                                 y_d_tomorrow            = n_d_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_d             = c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-1);
                                 exp_pc_d             = (y_d_tomorrow/c_d_tomorrow^(P.chic))*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow;
                                 exp_ee_d_jacob = exp_ee_d;

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_d_dax_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow*...
                                                (c_d_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_d_tomorrow) - P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_d_tomorrow);           

                                 % Derivatives with respect to b coefficients
                                 dpc_d_dbx_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_d_tomorrow,pi_tilde_d_tomorrow)*(pi_tilde_d_tomorrow-1)*pi_tilde_d_tomorrow) + y_d_tomorrow*(2*pi_tilde_d_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));                               


                                 exp_ee_d_int_temp         = exp_ee_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_d);
                                 exp_pc_d_int_temp         = exp_pc_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_d);

                                 dpc_d_dax_d_int = dpc_d_dax_d_int + pi^(-0.5)*G.e_weight(j)*dpc_d_dax_d;
                                 dpc_d_dbx_d_int = dpc_d_dbx_d_int + pi^(-0.5)*G.e_weight(j)*dpc_d_dbx_d;
                             else
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_d,C.T,C.P);

                                 %Deflationary Regime
                                 c_d_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_d,C.max,C.T,C.P);
                                 inf_d_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow_zlb     = inf_d_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow_zlb            = (c_d_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_d_tomorrow_zlb-1)^2));
                                 y_d_tomorrow_zlb            = n_d_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_d             = c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-1);
                                 exp_pc_zlb_d             = (y_d_tomorrow_zlb/c_d_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb;
                                 exp_ee_d_jacob = exp_ee_zlb_d;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dee_d_daxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today*P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_d_tomorrow_zlb^(-1);                     

                                 % Derivatives with respect to b coefficients
                                 dee_d_dbxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today*c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);   

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpc_d_daxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb*...
                                                (c_d_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_d_tomorrow_zlb) - P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_d_tomorrow_zlb);           

                                 % Derivatives with respect to b coefficients
                                 dpc_d_dbxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_d_tomorrow_zlb,pi_tilde_d_tomorrow_zlb)*(pi_tilde_d_tomorrow_zlb-1)*pi_tilde_d_tomorrow_zlb) + y_d_tomorrow_zlb*(2*pi_tilde_d_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                               

                                 exp_ee_d_int_temp         = exp_ee_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_d);
                                 exp_pc_d_int_temp         = exp_pc_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_d);

                                 dee_d_daxzlb_d_int = dee_d_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dee_d_daxzlb_d;
                                 dee_d_dbxzlb_d_int = dee_d_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dee_d_dbxzlb_d;

                                 dpc_d_daxzlb_d_int = dpc_d_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpc_d_daxzlb_d;
                                 dpc_d_dbxzlb_d_int = dpc_d_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpc_d_dbxzlb_d;

                             end              
                         end
                     end
                     exp_ee_d_int = cond_prob_sun*[exp_ee_s_int_temp;exp_ee_d_int_temp];
                     exp_pc_d_int = cond_prob_sun*[exp_pc_s_int_temp;exp_pc_d_int_temp];
                     
                     if r_d_tomorrow >= 1
                         % RHS derivative of EE (for Jacobian)
                         % Derivatives with respect to a coefficients
                         dee_d_dax_d = P.beta*del_today*(dr_dax(C.basis(:,i),inf_today,y_today,pi_tilde_today)*(cond_prob_sun*[exp_ee_s_jacob;exp_ee_d_jacob])...
                                     -cond_prob_sun(2)*r_today*P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_d_tomorrow^(-1));

                         % Derivatives with respect to b coefficients
                         dee_d_dbx_d = P.beta*del_today*(dr_dbx(C.basis(:,i),inf_today,y_today,c_today,pi_tilde_today)*(cond_prob_sun*[exp_ee_s_jacob;exp_ee_d_jacob])...
                                     - cond_prob_sun(2)*r_today*c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-2)*dpi_dbx(basis_interp));        

                         dee_d_dax_d_int = dee_d_dax_d_int + pi^(-0.5)*G.e_weight(j)*dee_d_dax_d;
                         dee_d_dbx_d_int = dee_d_dbx_d_int + pi^(-0.5)*G.e_weight(j)*dee_d_dbx_d;         
                     end
                end
                 % Residuals
                 res1_d(i) = c_today^(-P.chic) - P.beta*del_today*r_today*exp_ee_d_int;
                 res2_d(i) = (y_today/c_today^(P.chic))*(P.varphi*(pi_tilde_today-1)*pi_tilde_today-(1-P.theta)-P.theta*(1-P.tau)*w_today) - P.beta*del_today*exp_pc_d_int;
                 
                 % Jacobian
                 % Jacobian enteries relating to EE
                 lhs_jacob_ee_d(i,:) = [0*-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i)), 0*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))...
                                      -P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i)), 0*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))];
                 rhs_jacob_ee_d(i,:) = [dee_d_dax_s_int, dee_d_dbx_s_int, dee_d_daxzlb_s_int, dee_d_dbxzlb_s_int  ...
                                      dee_d_dax_d_int, dee_d_dbx_d_int, dee_d_daxzlb_d_int, dee_d_dbxzlb_d_int];

                 % Jacobian enteries relating to PC
                 lhs_jacob_pc_ax = (dy_dax(C.basis(:,i),pi_tilde_today)*c_today^(-P.chic) - P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i))*y_today)* ...
                                  (P.varphi*(pi_tilde_today - 1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
                 lhs_jacob_pc_ax = lhs_jacob_pc_ax + -P.theta*(1-P.tau)*dw_dax(C.basis(:,i),c_today,pi_tilde_today,n_today)* ...
                                  y_today/(c_today^(P.chic));

                 lhs_jacob_pc_bx = dy_dbx(C.basis(:,i),c_today,pi_tilde_today)/c_today^(P.chic)*(P.varphi*(pi_tilde_today-1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
                 lhs_jacob_pc_bx = lhs_jacob_pc_bx + y_today/(c_today^(P.chic))*(P.varphi*(2*pi_tilde_today*dpitilde_dbx(C.basis(:,i))-dpitilde_dbx(C.basis(:,i))) - P.theta*(1-P.tau)*dw_dbx(C.basis(:,i),c_today,pi_tilde_today,n_today));

                 lhs_jacob_pc_d(i,:)  = [0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx, 0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx...
                                       lhs_jacob_pc_ax, lhs_jacob_pc_bx, 0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx];
                 rhs_jacob_pc_d(i,:)  = [dpc_d_dax_s_int,dpc_d_dbx_s_int, dpc_d_daxzlb_s_int,dpc_d_dbxzlb_s_int...
                                       dpc_d_dax_d_int,dpc_d_dbx_d_int, dpc_d_daxzlb_d_int,dpc_d_dbxzlb_d_int];                
             else
                 c_today_zlb        = pf_c_zlb_d(i);
                 inf_today_zlb      = pf_inf_zlb_d(i);

                 pi_tilde_today_zlb     = inf_today_zlb/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
                 n_today_zlb           = (c_today_zlb/(1-(P.varphi/2)*(pi_tilde_today_zlb-1)^2));
                 y_today_zlb            = n_today_zlb;
                 w_today_zlb            = n_today_zlb^P.chin*c_today_zlb^P.chic;
                 r_today_zlb            = S.r_zlb;

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                 %For GH Integration
                exp_ee_s_int_zlb_temp = 0;
                exp_pc_s_int_zlb_temp = 0;
                exp_ee_d_int_zlb_temp = 0;
                exp_pc_d_int_zlb_temp = 0;
                exp_ee_d_int_zlb = 0;
                exp_pc_d_int_zlb = 0;

                deezlb_d_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_d_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_d_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_d_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_d_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_d_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_d_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                deezlb_d_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                dpczlb_d_dax_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_d_dax_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_d_daxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_d_daxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_d_dbx_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_d_dbx_d_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_d_dbxzlb_s_int    = zeros(1,O.n1+1); % (for Jacobian)
                dpczlb_d_dbxzlb_d_int    = zeros(1,O.n1+1); % (for Jacobian)

                for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_s,C.max,C.T,C.P);
                                 inf_s_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow     = inf_s_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow            = (c_s_tomorrow/(1-(P.varphi/2)*(pi_tilde_s_tomorrow-1)^2));
                                 y_s_tomorrow            = n_s_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_s             = c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-1);
                                 exp_pc_s             = (y_s_tomorrow/c_s_tomorrow^(P.chic))*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow;
                                 
                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_d_dax_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_s_tomorrow^(-1);

                                 % Derivatives with respect to b coefficients
                                 deezlb_d_dbx_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-2)*dpi_dbx(basis_interp);

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_d_dax_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow*...
                                                (c_s_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_s_tomorrow) - P.chic*c_s_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_s_tomorrow);

                                % Derivatives with respect to b coefficients
                                 dpczlb_d_dbx_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_s_tomorrow,pi_tilde_s_tomorrow)*(pi_tilde_s_tomorrow-1)*pi_tilde_s_tomorrow) + y_s_tomorrow*(2*pi_tilde_s_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));                    

                                 exp_ee_s_int_zlb_temp         = exp_ee_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_s);
                                 exp_pc_s_int_zlb_temp         = exp_pc_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_s);

                                 deezlb_d_dax_s_int = deezlb_d_dax_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_dax_s;
                                 deezlb_d_dbx_s_int = deezlb_d_dbx_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_dbx_s;

                                 dpczlb_d_dax_s_int = dpczlb_d_dax_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_dax_s;
                                 dpczlb_d_dbx_s_int = dpczlb_d_dbx_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_dbx_s;

                             else 
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_s,C.T,C.P);

                                 % Get Interpolated Values
                                 %Standard Regime
                                 c_s_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_s,C.max,C.T,C.P);
                                 inf_s_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_s,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_s_tomorrow_zlb     = inf_s_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_s_tomorrow_zlb            = (c_s_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_s_tomorrow_zlb-1)^2));
                                 y_s_tomorrow_zlb            = n_s_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_s             = c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-1);
                                 exp_pc_zlb_s             = (y_s_tomorrow_zlb/c_s_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_d_daxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_s_tomorrow_zlb^(-1);

                                 % Derivatives with respect to b coefficients
                                 deezlb_d_dbxzlb_s = -cond_prob_sun(1)*P.beta*del_today*r_today_zlb*c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_d_daxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb*...
                                                (c_s_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_s_tomorrow_zlb) - P.chic*c_s_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_s_tomorrow_zlb);

                                 % Derivatives with respect to b coefficients
                                 dpczlb_d_dbxzlb_s = cond_prob_sun(1)*P.beta*del_today*P.varphi/(c_s_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_s_tomorrow_zlb,pi_tilde_s_tomorrow_zlb)*(pi_tilde_s_tomorrow_zlb-1)*pi_tilde_s_tomorrow_zlb) + y_s_tomorrow_zlb*(2*pi_tilde_s_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                    


                                 exp_ee_s_int_zlb_temp         = exp_ee_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_s);
                                 exp_pc_s_int_zlb_temp         = exp_pc_s_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_s);

                                 deezlb_d_daxzlb_s_int = deezlb_d_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_daxzlb_s;
                                 deezlb_d_dbxzlb_s_int = deezlb_d_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_dbxzlb_s;

                                 dpczlb_d_daxzlb_s_int = dpczlb_d_daxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_daxzlb_s;
                                 dpczlb_d_dbxzlb_s_int = dpczlb_d_dbxzlb_s_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_dbxzlb_s;
                             end 
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Basis
                                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac_d,C.T,C.P);

                                 % Get Interpolated Values
                                 %Deflationary Regime
                                 c_d_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_d,C.max,C.T,C.P);
                                 inf_d_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow     = inf_d_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow            = (c_d_tomorrow/(1-(P.varphi/2)*(pi_tilde_d_tomorrow-1)^2));
                                 y_d_tomorrow            = n_d_tomorrow;

                                 % Compute GH Integration
                                 exp_ee_d             = c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-1);
                                 exp_pc_d             = (y_d_tomorrow/c_d_tomorrow^(P.chic))*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_d_dax_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_d_tomorrow^(-1);                     

                                 % Derivatives with respect to b coefficients
                                 deezlb_d_dbx_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-2)*dpi_dbx(basis_interp);   

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_d_dax_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow*...
                                                (c_d_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_d_tomorrow) - P.chic*c_d_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_d_tomorrow);           

                                 % Derivatives with respect to b coefficients
                                 dpczlb_d_dbx_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow^P.chic)*...
                                                ((dy_dbx(basis_interp,c_d_tomorrow,pi_tilde_d_tomorrow)*(pi_tilde_d_tomorrow-1)*pi_tilde_d_tomorrow) + y_d_tomorrow*(2*pi_tilde_d_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));                               

                                 exp_ee_d_int_zlb_temp         = exp_ee_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_d);
                                 exp_pc_d_int_zlb_temp         = exp_pc_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_d);

                                 deezlb_d_dax_d_int = deezlb_d_dax_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_dax_d;
                                 deezlb_d_dbx_d_int = deezlb_d_dbx_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_dbx_d;

                                 dpczlb_d_dax_d_int = dpczlb_d_dax_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_dax_d;
                                 dpczlb_d_dbx_d_int = dpczlb_d_dbx_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_dbx_d;
                             else
                                 % Get Interpolated Basis
                                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb_d,C.T,C.P);

                                 %Deflationary Regime
                                 c_d_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac_zlb_d,C.max,C.T,C.P);
                                 inf_d_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf_zlb_d,C.max,C.T,C.P);

                                 % Build out other PF Values
                                 pi_tilde_d_tomorrow_zlb     = inf_d_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                                 n_d_tomorrow_zlb            = (c_d_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_d_tomorrow_zlb-1)^2));
                                 y_d_tomorrow_zlb            = n_d_tomorrow_zlb;

                                 % Compute GH Integration
                                 exp_ee_zlb_d             = c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-1);
                                 exp_pc_zlb_d             = (y_d_tomorrow_zlb/c_d_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb;

                                 % RHS derivative of EE (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 deezlb_d_daxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_d_tomorrow_zlb^(-1);                     

                                 % Derivatives with respect to b coefficients
                                 deezlb_d_dbxzlb_d = -cond_prob_sun(2)*P.beta*del_today*r_today_zlb*c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb);   

                                 % RHS derivative of PC (for Jacobian)
                                 % Derivatives with respect to a coefficients
                                 dpczlb_d_daxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb*...
                                                (c_d_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_d_tomorrow_zlb) - P.chic*c_d_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_d_tomorrow_zlb);           

                                 % Derivatives with respect to b coefficients
                                 dpczlb_d_dbxzlb_d = cond_prob_sun(2)*P.beta*del_today*P.varphi/(c_d_tomorrow_zlb^P.chic)*...
                                                ((dy_dbx(basis_interp_zlb,c_d_tomorrow_zlb,pi_tilde_d_tomorrow_zlb)*(pi_tilde_d_tomorrow_zlb-1)*pi_tilde_d_tomorrow_zlb) + y_d_tomorrow_zlb*(2*pi_tilde_d_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));                               

                                 exp_ee_d_int_zlb_temp         = exp_ee_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_d);
                                 exp_pc_d_int_zlb_temp         = exp_pc_d_int_zlb_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_d);

                                 deezlb_d_daxzlb_d_int = deezlb_d_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_daxzlb_d;
                                 deezlb_d_dbxzlb_d_int = deezlb_d_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*deezlb_d_dbxzlb_d;

                                 dpczlb_d_daxzlb_d_int = dpczlb_d_daxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_daxzlb_d;
                                 dpczlb_d_dbxzlb_d_int = dpczlb_d_dbxzlb_d_int + pi^(-0.5)*G.e_weight(j)*dpczlb_d_dbxzlb_d; 

                             end              
                         end
                     end
                     exp_ee_d_int_zlb = cond_prob_sun*[exp_ee_s_int_zlb_temp;exp_ee_d_int_zlb_temp];
                     exp_pc_d_int_zlb = cond_prob_sun*[exp_pc_s_int_zlb_temp;exp_pc_d_int_zlb_temp];         
                 end

                 % Residuals
                 res1_zlb_d(i) = c_today_zlb^(-P.chic) - P.beta*del_today*r_today_zlb*exp_ee_d_int_zlb;
                 res2_zlb_d(i) = (y_today_zlb/c_today_zlb^(P.chic))*(P.varphi*(pi_tilde_today_zlb-1)*pi_tilde_today_zlb-(1-P.theta)-P.theta*(1-P.tau)*w_today_zlb) - P.beta*del_today*exp_pc_d_int_zlb;
                 
                 % Jacobian
                 % Jacobian enteries relating to EE
                 lhs_jacob_ee_zlb_d(i,:) = [0*-P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i)), 0*-P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))...
                                    0*-P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i)), -P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i))];
                 rhs_jacob_ee_zlb_d(i,:) = [deezlb_d_dax_s_int, deezlb_d_dbx_s_int, deezlb_d_daxzlb_s_int, deezlb_d_dbxzlb_s_int ...
                                      deezlb_d_dax_d_int, deezlb_d_dbx_d_int, deezlb_d_daxzlb_d_int, deezlb_d_dbxzlb_d_int];

                 % Jacobian enteries relating to PC
                 lhs_jacob_pc_ax = (dy_dax(C.basis(:,i),pi_tilde_today_zlb)*c_today_zlb^(-P.chic) - P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i))*y_today_zlb)* ...
                                  (P.varphi*(pi_tilde_today_zlb - 1)*pi_tilde_today_zlb - (1-P.theta) - P.theta*(1-P.tau)*w_today_zlb);
                 lhs_jacob_pc_ax = lhs_jacob_pc_ax + -P.theta*(1-P.tau)*dw_dax(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb,n_today_zlb)* ...
                                  y_today_zlb/(c_today_zlb^(P.chic));

                 lhs_jacob_pc_bx = dy_dbx(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb)/c_today_zlb^(P.chic)*(P.varphi*(pi_tilde_today_zlb-1)*pi_tilde_today_zlb - (1-P.theta) - P.theta*(1-P.tau)*w_today_zlb);
                 lhs_jacob_pc_bx = lhs_jacob_pc_bx + y_today_zlb/(c_today_zlb^(P.chic))*(P.varphi*(2*pi_tilde_today_zlb*dpitilde_dbx(C.basis(:,i))-dpitilde_dbx(C.basis(:,i))) - P.theta*(1-P.tau)*dw_dbx(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb,n_today_zlb));

                 lhs_jacob_pc_zlb_d(i,:)  = [0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx, 0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx...
                                       0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx, lhs_jacob_pc_ax, lhs_jacob_pc_bx];
                 rhs_jacob_pc_zlb_d(i,:)  = [dpczlb_d_dax_s_int,dpczlb_d_dbx_s_int, dpczlb_d_daxzlb_s_int,dpczlb_d_dbxzlb_s_int...
                                       dpczlb_d_dax_d_int,dpczlb_d_dbx_d_int, dpczlb_d_daxzlb_d_int,dpczlb_d_dbxzlb_d_int];
             end
         end
     end     
 end
 
 R = [res1_s,res2_s,res1_zlb_s,res2_zlb_s,res1_d,res2_d,res1_zlb_d,res2_zlb_d];
 
 J_nzlb_s = [lhs_jacob_ee_s-rhs_jacob_ee_s;lhs_jacob_pc_s-rhs_jacob_pc_s];
 J_zlb_s = [lhs_jacob_ee_zlb_s-rhs_jacob_ee_zlb_s;lhs_jacob_pc_zlb_s-rhs_jacob_pc_zlb_s]; 
 J_nzlb_d =  [lhs_jacob_ee_d-rhs_jacob_ee_d;lhs_jacob_pc_d-rhs_jacob_pc_d];
 J_zlb_d =  [lhs_jacob_ee_zlb_d-rhs_jacob_ee_zlb_d;lhs_jacob_pc_zlb_d-rhs_jacob_pc_zlb_d];
 
 J = [J_nzlb_s;J_zlb_s;J_nzlb_d;J_zlb_d];