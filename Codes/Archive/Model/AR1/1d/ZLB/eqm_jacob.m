function [R,J] = eqm_jacob(weights,P,S,G,O,C)

% [R,J] = eqm_jacob(weights,P,S,G,O,C)
%   Outputs residuals and Analytical Jacobian of the equilibrium system of 
%   equations for least squares minimization/Chebyshev interpolation method.
% Inputs:
%   weights :   Parameter Weights for Cheb Poly (Matrix)
%   P       :   Structure of parameters
%   S       :   Structure of steady state values
%   G       :   Structure of grids
%   O       :   Structure of options
%   C       :   Structure of Chebyshev polynomial parameters
% Output:
%   R       :   Residuals
%   J       :   Jacobian

global pi_yesterday

% Allocate Space for Residuals

res1        = zeros(G.nodes,1);
res2        = zeros(G.nodes,1);
res1_zlb        = zeros(G.nodes,1);
res2_zlb        = zeros(G.nodes,1);

[~,col] = size(weights);

lhs_jacob_ee = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_ee = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_pc = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_pc = zeros(G.nodes,col*(O.n1+1));


lhs_jacob_ee_zlb = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_ee_zlb = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_pc_zlb = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_pc_zlb = zeros(G.nodes,col*(O.n1+1));


% Map Coefficients 
C_Ac    = weights(1:O.n1+1,1);
C_Ainf  = weights(1:O.n1+1,2);
C_Ac_zlb    = weights(1:O.n1+1,3);
C_Ainf_zlb  = weights(1:O.n1+1,4);

% Get Partial Derivatives for Analytical Jacobian
get_partial_derivs;

% Get Policy Functions
pf_c = allcheb111(O.delbound,G.del_grid,C_Ac,C.T,C.P);    
pf_inf = allcheb111(O.delbound,G.del_grid,C_Ainf,C.T,C.P); 

pf_pitilde  = pf_inf/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
pf_n        = (pf_c./(1-(P.varphi/2).*(pf_pitilde-1).^2));
pf_y        = pf_n;
pf_w        = pf_n.^P.chin.*pf_c.^P.chic;
pf_r        = S.pi_targ/P.beta*((pf_inf./S.pi_targ).^(P.phi_pi).*(pf_y./S.y).^(P.phi_y));


pf_c_zlb = allcheb111(O.delbound,G.del_grid,C_Ac_zlb,C.T,C.P);    
pf_inf_zlb = allcheb111(O.delbound,G.del_grid,C_Ainf_zlb,C.T,C.P); 
pf_r_zlb = ones(size(pf_r));
 
 for i = 1:G.nodes
     del_today         = G.del_gr(i);
     
     c_temp            = pf_c(i);
     inf_temp          = pf_inf(i);

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
         
         % Get r_today's cheb weights
         
         C_Ar = chebweights11(O.n1,pf_r,C.T,C.P,C.X);

         %For GH Integration
         exp_ee_int = 0;
         exp_pc_int = 0;
         
         dee_dax_int    = zeros(1,O.n1+1); % (for Jacobian)
         dee_daxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)
         dee_dbx_int    = zeros(1,O.n1+1); % (for Jacobian)
         dee_dbxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)
         
         dpc_dax_int    = zeros(1,O.n1+1); % (for Jacobian)
         dpc_daxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)
         dpc_dbx_int    = zeros(1,O.n1+1); % (for Jacobian)
         dpc_dbxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)

         for j = 1:length(G.e_weight)
             r_tomorrow = allcheb111(O.delbound,del_tomorrow(j),C_Ar,C.T,C.P);
             
             if r_tomorrow >= 1
                 % Get Interpolated Basis
                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac,C.T,C.P);

                 % Get Interpolated Values
                 c_tomorrow         = allcheb111(O.delbound,del_tomorrow(j),C_Ac,C.T,C.P);
                 inf_tomorrow       = allcheb111(O.delbound,del_tomorrow(j),C_Ainf,C.T,C.P);

                 % Build out other PF Values
                 pi_tilde_tomorrow     = inf_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                 n_tomorrow         = (c_tomorrow/(1-(P.varphi/2)*(pi_tilde_tomorrow-1)^2));
                 y_tomorrow         = n_tomorrow;

                 exp_ee = c_tomorrow^(-P.chic)*inf_tomorrow^(-1);
                 exp_pc = (y_tomorrow/c_tomorrow^(P.chic))*P.varphi*(pi_tilde_tomorrow - 1)*pi_tilde_tomorrow;

                 % RHS derivative of EE (for Jacobian)
                 % Derivatives with respect to a coefficients
                 dee_dax = P.beta*del_today*(dr_dax(C.basis(:,i),inf_today,y_today,pi_tilde_today)*(c_tomorrow^(-P.chic)*inf_tomorrow^(-1))...
                             - r_today*P.chic*c_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_tomorrow^(-1));             

                 % Derivatives with respect to b coefficients
                 dee_dbx = P.beta*del_today*(dr_dbx(C.basis(:,i),inf_today,y_today,c_today,pi_tilde_today)*(c_tomorrow^(-P.chic)*inf_tomorrow^(-1))...
                             - r_today*c_tomorrow^(-P.chic)*inf_tomorrow^(-2)*dpi_dbx(basis_interp));

                 % RHS derivative of PC (for Jacobian)
                 % Derivatives with respect to a coefficients
                 dpc_dax = P.beta*del_today*P.varphi*(pi_tilde_tomorrow - 1)*pi_tilde_tomorrow*...
                                (c_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_tomorrow) - P.chic*c_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_tomorrow);

                 % Derivatives with respect to b coefficients
                 dpc_dbx = P.beta*del_today*P.varphi/(c_tomorrow^P.chic)*...
                                ((dy_dbx(basis_interp,c_tomorrow,pi_tilde_tomorrow)*(pi_tilde_tomorrow-1)*pi_tilde_tomorrow) + y_tomorrow*(2*pi_tilde_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));

                 
                 % Compute GH Integration
                 exp_ee_int = exp_ee_int + pi^(-0.5)*G.e_weight(j)*(exp_ee);
                 exp_pc_int = exp_pc_int + pi^(-0.5)*G.e_weight(j)*(exp_pc);

                 dee_dax_int = dee_dax_int + pi^(-0.5)*G.e_weight(j)*dee_dax;
                 dee_dbx_int = dee_dbx_int + pi^(-0.5)*G.e_weight(j)*dee_dbx;
                 
                 dpc_dax_int = dpc_dax_int + pi^(-0.5)*G.e_weight(j)*dpc_dax;
                 dpc_dbx_int = dpc_dbx_int + pi^(-0.5)*G.e_weight(j)*dpc_dbx;

             else 
                 % Get Interpolated Basis
                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb,C.T,C.P);
                 
                 % Get Interpolated Values
                 c_tomorrow_zlb     = allcheb111(O.delbound,del_tomorrow(j),C_Ac_zlb ,C.T,C.P);
                 inf_tomorrow_zlb   = allcheb111(O.delbound,del_tomorrow(j),C_Ainf_zlb ,C.T,C.P);
                 % Build out other PF Values
                 pi_tilde_tomorrow_zlb     = inf_tomorrow_zlb/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
                 n_tomorrow_zlb         = (c_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_tomorrow_zlb-1)^2));
                 y_tomorrow_zlb         = n_tomorrow_zlb;

                 exp_ee_zlb = c_tomorrow_zlb^(-P.chic)*inf_tomorrow_zlb^(-1);
                 exp_pc_zlb = (y_tomorrow_zlb/c_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_tomorrow_zlb - 1)*pi_tilde_tomorrow_zlb;

                 % RHS derivative of EE (for Jacobian)
                 % Derivatives with respect to a coefficients
                 dee_daxzlb = -P.beta*del_today*(r_today*P.chic*c_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_tomorrow_zlb^(-1));             

                 % Derivatives with respect to b coefficients
                 dee_dbxzlb = -P.beta*del_today*(r_today*c_tomorrow_zlb^(-P.chic)*inf_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb));

                 % RHS derivative of PC (for Jacobian)
                 % Derivatives with respect to a coefficients
                 dpc_daxzlb = P.beta*del_today*P.varphi*(pi_tilde_tomorrow_zlb - 1)*pi_tilde_tomorrow_zlb*...
                                (c_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_tomorrow_zlb) - P.chic*c_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_tomorrow_zlb);

                 % Derivatives with respect to b coefficients
                 dpc_dbxzlb = P.beta*del_today*P.varphi/(c_tomorrow_zlb^P.chic)*...
                                ((dy_dbx(basis_interp_zlb,c_tomorrow_zlb,pi_tilde_tomorrow_zlb)*(pi_tilde_tomorrow_zlb-1)*pi_tilde_tomorrow_zlb) + y_tomorrow_zlb*(2*pi_tilde_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));

                 % Compute GH Integration
                 exp_ee_int = exp_ee_int + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb);
                 exp_pc_int = exp_pc_int + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb);

                 dee_daxzlb_int = dee_daxzlb_int + pi^(-0.5)*G.e_weight(j)*dee_daxzlb;
                 dee_dbxzlb_int = dee_dbxzlb_int + pi^(-0.5)*G.e_weight(j)*dee_dbxzlb;

                 dpc_daxzlb_int = dpc_daxzlb_int + pi^(-0.5)*G.e_weight(j)*dpc_daxzlb;
                 dpc_dbxzlb_int = dpc_dbxzlb_int + pi^(-0.5)*G.e_weight(j)*dpc_dbxzlb;
             end            
         end
         % Residuals
         res1(i) = c_today^(-P.chic) - P.beta*del_today*r_today*exp_ee_int;
         res2(i) = (y_today/c_today^(P.chic))*(P.varphi*(pi_tilde_today-1)*pi_tilde_today-(1-P.theta)-P.theta*(1-P.tau)*w_today) - P.beta*del_today*exp_pc_int;
         
         % Jacobian
         % Jacobian enteries relating to EE
         lhs_jacob_ee(i,:) = [-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i))...
                            0*-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))];
         rhs_jacob_ee(i,:) = [dee_dax_int, dee_dbx_int ...
                            dee_daxzlb_int, dee_dbxzlb_int];

         % Jacobian enteries relating to PC
         lhs_jacob_pc_ax = (dy_dax(C.basis(:,i),pi_tilde_today)*c_today^(-P.chic) - P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i))*y_today)* ...
                          (P.varphi*(pi_tilde_today - 1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
         lhs_jacob_pc_ax = lhs_jacob_pc_ax + -P.theta*(1-P.tau)*dw_dax(C.basis(:,i),c_today,pi_tilde_today,n_today)* ...
                          y_today/(c_today^(P.chic));

         lhs_jacob_pc_bx = dy_dbx(C.basis(:,i),c_today,pi_tilde_today)/c_today^(P.chic)*(P.varphi*(pi_tilde_today-1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
         lhs_jacob_pc_bx = lhs_jacob_pc_bx + y_today/(c_today^(P.chic))*(P.varphi*(2*pi_tilde_today*dpitilde_dbx(C.basis(:,i))-dpitilde_dbx(C.basis(:,i))) - P.theta*(1-P.tau)*dw_dbx(C.basis(:,i),c_today,pi_tilde_today,n_today));

         lhs_jacob_pc(i,:)  = [lhs_jacob_pc_ax, lhs_jacob_pc_bx...
                               0*lhs_jacob_pc_ax, 0*lhs_jacob_pc_bx];
         rhs_jacob_pc(i,:)  = [dpc_dax_int,dpc_dbx_int...
                               dpc_daxzlb_int,dpc_dbxzlb_int];                       
     else
         c_today_zlb        = pf_c_zlb(i);
         inf_today_zlb      = pf_inf_zlb(i);

         pi_tilde_today_zlb     = inf_today_zlb/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
         n_today_zlb           = (c_today_zlb/(1-(P.varphi/2)*(pi_tilde_today_zlb-1)^2));
         y_today_zlb            = n_today_zlb;
         w_today_zlb            = n_today_zlb^P.chin*c_today_zlb^P.chic;
         r_today_zlb            = S.r_zlb;
              
         % Get Tomorrow's Shock 
         del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;
         
         % Get r_today's cheb weights
         
         C_Ar = chebweights11(O.n1,pf_r,C.T,C.P,C.X);

         %For GH Integration
         exp_ee_int_zlb = 0;
         exp_pc_int_zlb = 0;
         
         deezlb_dax_int    = zeros(1,O.n1+1); % (for Jacobian)
         deezlb_daxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)
         deezlb_dbx_int    = zeros(1,O.n1+1); % (for Jacobian)
         deezlb_dbxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)
         
         dpczlb_dax_int    = zeros(1,O.n1+1); % (for Jacobian)
         dpczlb_daxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)
         dpczlb_dbx_int    = zeros(1,O.n1+1); % (for Jacobian)
         dpczlb_dbxzlb_int    = zeros(1,O.n1+1); % (for Jacobian)

         for j = 1:length(G.e_weight)
             r_tomorrow = allcheb111(O.delbound,del_tomorrow(j),C_Ar,C.T,C.P);
             
             if r_tomorrow >= 1
                 % Get Interpolated Basis
                 basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac,C.T,C.P);

                 % Get Interpolated Values
                 c_tomorrow         = allcheb111(O.delbound,del_tomorrow(j),C_Ac,C.T,C.P);
                 inf_tomorrow       = allcheb111(O.delbound,del_tomorrow(j),C_Ainf,C.T,C.P);
                 
                 % Build out other PF Values
                 pi_tilde_tomorrow     = inf_tomorrow/(S.pi_targ^P.iota*inf_today_zlb^(1-P.iota))^P.alpha;
                 n_tomorrow         = (c_tomorrow/(1-(P.varphi/2)*(pi_tilde_tomorrow-1)^2));
                 y_tomorrow         = n_tomorrow;

                 exp_ee = c_tomorrow^(-P.chic)*inf_tomorrow^(-1);
                 exp_pc = (y_tomorrow/c_tomorrow^(P.chic))*P.varphi*(pi_tilde_tomorrow - 1)*pi_tilde_tomorrow;

                 % RHS derivative of EE (for Jacobian)
                 % Derivatives with respect to a coefficients
                 deezlb_dax = -P.beta*del_today*r_today_zlb*P.chic*c_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_tomorrow^(-1);             

                 % Derivatives with respect to b coefficients
                 deezlb_dbx = -P.beta*del_today*r_today_zlb*c_tomorrow^(-P.chic)*inf_tomorrow^(-2)*dpi_dbx(basis_interp);

                 % RHS derivative of PC (for Jacobian)
                 % Derivatives with respect to a coefficients
                 dpczlb_dax = P.beta*del_today*P.varphi*(pi_tilde_tomorrow - 1)*pi_tilde_tomorrow*...
                                (c_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_tomorrow) - P.chic*c_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_tomorrow);

                 % Derivatives with respect to b coefficients
                 dpczlb_dbx = P.beta*del_today*P.varphi/(c_tomorrow^P.chic)*...
                                ((dy_dbx(basis_interp,c_tomorrow,pi_tilde_tomorrow)*(pi_tilde_tomorrow-1)*pi_tilde_tomorrow) + y_tomorrow*(2*pi_tilde_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));

                 % Compute GH Integration
                 exp_ee_int_zlb = exp_ee_int_zlb + pi^(-0.5)*G.e_weight(j)*(exp_ee);
                 exp_pc_int_zlb = exp_pc_int_zlb + pi^(-0.5)*G.e_weight(j)*(exp_pc);

                 deezlb_dax_int = deezlb_dax_int + pi^(-0.5)*G.e_weight(j)*deezlb_dax;
                 deezlb_dbx_int = deezlb_dbx_int + pi^(-0.5)*G.e_weight(j)*deezlb_dbx;

                 dpczlb_dax_int = dpczlb_dax_int + pi^(-0.5)*G.e_weight(j)*dpczlb_dax;
                 dpczlb_dbx_int = dpczlb_dbx_int + pi^(-0.5)*G.e_weight(j)*dpczlb_dbx;
             else 
                 % Get Interpolated Basis
                 basis_interp_zlb = basisinterp(O.delbound,del_tomorrow(j),C_Ac_zlb,C.T,C.P);
                 
                 % Get Interpolated Values
                 c_tomorrow_zlb     = allcheb111(O.delbound,del_tomorrow(j),C_Ac_zlb ,C.T,C.P);
                 inf_tomorrow_zlb   = allcheb111(O.delbound,del_tomorrow(j),C_Ainf_zlb ,C.T,C.P);

                 % Build out other PF Values
                 pi_tilde_tomorrow_zlb     = inf_tomorrow_zlb/(S.pi_targ^P.iota*inf_today_zlb^(1-P.iota))^P.alpha;
                 n_tomorrow_zlb         = (c_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_tomorrow_zlb-1)^2));
                 y_tomorrow_zlb         = n_tomorrow_zlb;

                 exp_ee_zlb = c_tomorrow_zlb^(-P.chic)*inf_tomorrow_zlb^(-1);
                 exp_pc_zlb = (y_tomorrow_zlb/c_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_tomorrow_zlb - 1)*pi_tilde_tomorrow_zlb;

                 % RHS derivative of EE (for Jacobian)
                 % Derivatives with respect to a coefficients
                 deezlb_daxzlb = P.beta*del_today*(drzlb_dax(C.basis(:,i),inf_today_zlb,y_today_zlb,pi_tilde_today_zlb)*(c_tomorrow_zlb^(-P.chic)*inf_tomorrow_zlb^(-1))...
                             - r_today_zlb*P.chic*c_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*inf_tomorrow_zlb^(-1));

                 % Derivatives with respect to b coefficients
                 deezlb_dbxzlb = P.beta*del_today*(drzlb_dbx(C.basis(:,i),inf_today_zlb,y_today_zlb,c_today_zlb,pi_tilde_today_zlb)*(c_tomorrow_zlb^(-P.chic)*inf_tomorrow_zlb^(-1))...
                             - r_today_zlb*c_tomorrow_zlb^(-P.chic)*inf_tomorrow_zlb^(-2)*dpi_dbx(basis_interp_zlb));

                 % RHS derivative of PC (for Jacobian)
                 % Derivatives with respect to a coefficients
                 dpczlb_daxzlb = P.beta*del_today*P.varphi*(pi_tilde_tomorrow_zlb - 1)*pi_tilde_tomorrow_zlb*...
                                (c_tomorrow_zlb^(-P.chic)*dy_dax(basis_interp_zlb,pi_tilde_tomorrow_zlb) - P.chic*c_tomorrow_zlb^(-P.chic-1)*dc_dax(basis_interp_zlb)*y_tomorrow_zlb);

                 % Derivatives with respect to b coefficients
                 dpczlb_dbxzlb = P.beta*del_today*P.varphi/(c_tomorrow_zlb^P.chic)*...
                                ((dy_dbx(basis_interp_zlb,c_tomorrow_zlb,pi_tilde_tomorrow_zlb)*(pi_tilde_tomorrow_zlb-1)*pi_tilde_tomorrow_zlb) + y_tomorrow_zlb*(2*pi_tilde_tomorrow_zlb*dpitilde_dbx(basis_interp_zlb) - dpitilde_dbx(basis_interp_zlb)));

                 % Compute GH Integration
                 exp_ee_int_zlb = exp_ee_int_zlb + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb);
                 exp_pc_int_zlb = exp_pc_int_zlb + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb);

                 deezlb_daxzlb_int = deezlb_daxzlb_int + pi^(-0.5)*G.e_weight(j)*deezlb_daxzlb;
                 deezlb_dbxzlb_int = deezlb_dbxzlb_int + pi^(-0.5)*G.e_weight(j)*deezlb_dbxzlb;

                 dpczlb_daxzlb_int = dpczlb_daxzlb_int + pi^(-0.5)*G.e_weight(j)*dpczlb_daxzlb;
                 dpczlb_dbxzlb_int = dpczlb_dbxzlb_int + pi^(-0.5)*G.e_weight(j)*dpczlb_dbxzlb;
             end            
         end
        % Residuals
         res1_zlb(i) = c_today_zlb^(-P.chic) - P.beta*del_today*r_today_zlb*exp_ee_int_zlb;
         res2_zlb(i) = (y_today_zlb/c_today_zlb^(P.chic))*(P.varphi*(pi_tilde_today_zlb-1)*pi_tilde_today_zlb-(1-P.theta)-P.theta*(1-P.tau)*w_today_zlb) - P.beta*del_today*exp_pc_int_zlb;

         % Jacobian
         % Jacobian enteries relating to EE
         lhs_jacob_ee_zlb(i,:) = [0*-P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), 0*dc_dbx(C.basis(:,i))...
                              -P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i))];
         rhs_jacob_ee_zlb(i,:) = [deezlb_dax_int, deezlb_dbx_int ...
                            deezlb_daxzlb_int, deezlb_dbxzlb_int];

         % Jacobian enteries relating to PC
         lhs_jacob_pc_ax_zlb = (dy_dax(C.basis(:,i),pi_tilde_today_zlb)*c_today_zlb^(-P.chic) - P.chic*c_today_zlb^(-P.chic-1)*dc_dax(C.basis(:,i))*y_today_zlb)* ...
                          (P.varphi*(pi_tilde_today_zlb - 1)*pi_tilde_today_zlb - (1-P.theta) - P.theta*(1-P.tau)*w_today_zlb);
         lhs_jacob_pc_ax_zlb = lhs_jacob_pc_ax_zlb + -P.theta*(1-P.tau)*dw_dax(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb,n_today_zlb)* ...
                          y_today_zlb/(c_today_zlb^(P.chic));

         lhs_jacob_pc_bx_zlb = dy_dbx(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb)/c_today_zlb^(P.chic)*(P.varphi*(pi_tilde_today_zlb-1)*pi_tilde_today_zlb - (1-P.theta) - P.theta*(1-P.tau)*w_today_zlb);
         lhs_jacob_pc_bx_zlb = lhs_jacob_pc_bx_zlb + y_today_zlb/(c_today_zlb^(P.chic))*(P.varphi*(2*pi_tilde_today_zlb*dpitilde_dbx(C.basis(:,i))-dpitilde_dbx(C.basis(:,i))) - P.theta*(1-P.tau)*dw_dbx(C.basis(:,i),c_today_zlb,pi_tilde_today_zlb,n_today_zlb));

         lhs_jacob_pc_zlb(i,:)  = [0*lhs_jacob_pc_ax_zlb, 0*lhs_jacob_pc_bx_zlb...
                               lhs_jacob_pc_ax_zlb, lhs_jacob_pc_bx_zlb];
         rhs_jacob_pc_zlb(i,:)  = [dpczlb_dax_int,dpczlb_dbx_int...
                               dpczlb_daxzlb_int,dpczlb_dbxzlb_int];
     end
 end
 
 R = [res1,res2,res1_zlb,res2_zlb];
 J_NZLB = [lhs_jacob_ee-rhs_jacob_ee;lhs_jacob_pc-rhs_jacob_pc];
 J_ZLB = [lhs_jacob_ee_zlb-rhs_jacob_ee_zlb;lhs_jacob_pc_zlb-rhs_jacob_pc_zlb];
 
 J = [J_NZLB;J_ZLB];
