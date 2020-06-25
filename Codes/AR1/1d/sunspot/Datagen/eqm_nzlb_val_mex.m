function [R,J] = eqm_nzlb_val_mex(weights,P,S,G,O,C)

% R = eqm(x,state,pf,P,S,G,tbind)
%   Outputs residuals of the equilibrium system of equations for time
%   iteration/linear interpolation method.
% Inputs:
%   weights :   Parameter Weights for Cheb Poly (Matrix)
%   P       :   Structure of parameters
%   S       :   Structure of steady state values
%   G       :   Structure of grids
%   O       :   Structure of options
%   C       :   Structure of Chebyshev polynomial parameters
% Output:
%   R       :   Residuals

% weights

% -------------------------------------------------------------------------
% House Keeping
% -------------------------------------------------------------------------
global pi_yesterday

% Allocate Space for Residuals
res1 = zeros(G.nodes,1);
res2 = zeros(G.nodes,1);
res3 = zeros(G.nodes,1);

[~,col] = size(weights);

lhs_jacob_ee = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_ee = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_pc = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_pc = zeros(G.nodes,col*(O.n1+1));
lhs_jacob_v  = zeros(G.nodes,col*(O.n1+1));
rhs_jacob_v  = zeros(G.nodes,col*(O.n1+1));

% Map Coefficients 
C_Ac    = weights(:,1);
C_Ainf  = weights(:,2);
C_Av    = weights(:,3);

% Get Partial Derivatives for Analytical Jacobian
get_partial_derivs;

% Get Policy Functions
pf_c = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Ac,C.max,C.T,C.P);     
pf_inf = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Ainf,C.max,C.T,C.P);  
pf_v = Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Av,C.max,C.T,C.P);     

% -------------------------------------------------------------------------
% Main Code
% -------------------------------------------------------------------------

for i = 1:G.nodes
     del_today          = G.del_grid(i);
     c_today            = pf_c(i);
     inf_today          = pf_inf(i);
     v_today            = pf_v(i);

     % Build out other PFs rom current values
     pi_tilde_today     = inf_today/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
     n_today            = (c_today/(1-(P.varphi/2)*(pi_tilde_today-1)^2));
     y_today            = n_today;
     w_today            = n_today^P.chin*c_today^P.chic;
     r_today            = S.pi_targ/P.beta*((inf_today/S.pi_targ)^(P.phi_pi)*(y_today/S.y)^(P.phi_y));

     % Get Tomorrow's Shock 
     del_tomorrow       = P.rho*(del_today - 1) + 1 + G.e_nodes;

     %For GH Integration
     exp_ee_int = 0;
     exp_pc_int = 0;
     exp_v_int  = 0;
     
     dee_dax_int    = zeros(1,O.n1+1); % (for Jacobian)
     dee_dbx_int    = zeros(1,O.n1+1); % (for Jacobian)
     dee_dcx_int    = zeros(1,O.n1+1); % (for Jacobian)
     
     dpc_dax_int    = zeros(1,O.n1+1); % (for Jacobian)
     dpc_dbx_int    = zeros(1,O.n1+1); % (for Jacobian)
     dpc_dcx_int    = zeros(1,O.n1+1); % (for Jacobian)
     
     dv_dax_int     = zeros(1,O.n1+1); % (for Jacobian)
     dv_dbx_int     = zeros(1,O.n1+1); % (for Jacobian)
     dv_dcx_int     = zeros(1,O.n1+1); % (for Jacobian)


     for j = 1:length(G.e_weight)
         % Get Interpolated Basis
         basis_interp = basisinterp(O.delbound,del_tomorrow(j),C_Ac,C.T,C.P); 

         % Get Interpolated Values
         c_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ac,C.max,C.T,C.P); 
         inf_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ainf,C.max,C.T,C.P); 
         v_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av,C.max,C.T,C.P); 
         
         % Build out other PF Values
         pi_tilde_tomorrow     = inf_tomorrow/(S.pi_targ^P.iota*inf_today^(1-P.iota))^P.alpha;
         n_tomorrow         = (c_tomorrow/(1-(P.varphi/2)*(pi_tilde_tomorrow-1)^2));
         y_tomorrow         = n_tomorrow;
         
         exp_ee = c_tomorrow^(-P.chic)*inf_tomorrow^(-1);
         exp_pc = (y_tomorrow/c_tomorrow^(P.chic))*P.varphi*(pi_tilde_tomorrow - 1)*pi_tilde_tomorrow;
         exp_v = v_tomorrow;

         % RHS derivative of EE (for Jacobian)
         % Derivatives with respect to a coefficients
         dee_dax = P.beta*del_today*(dr_dax(C.basis(:,i),inf_today,y_today,pi_tilde_today)*(c_tomorrow^(-P.chic)*inf_tomorrow^(-1))...
                     - r_today*P.chic*c_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*inf_tomorrow^(-1));             

         % Derivatives with respect to b coefficients
         dee_dbx = P.beta*del_today*(dr_dbx(C.basis(:,i),inf_today,y_today,c_today,pi_tilde_today)*(c_tomorrow^(-P.chic)*inf_tomorrow^(-1))...
                     - r_today*c_tomorrow^(-P.chic)*inf_tomorrow^(-2)*dpi_dbx(basis_interp));

         % Derivatives with respect to c and c_d coefficients
         dee_dcx = 0*basis_interp';

         % RHS derivative of PC (for Jacobian)
         % Derivatives with respect to a coefficients
         dpc_dax = P.beta*del_today*P.varphi*(pi_tilde_tomorrow - 1)*pi_tilde_tomorrow*...
                        (c_tomorrow^(-P.chic)*dy_dax(basis_interp,pi_tilde_tomorrow) - P.chic*c_tomorrow^(-P.chic-1)*dc_dax(basis_interp)*y_tomorrow);
                    
         % Derivatives with respect to b coefficients
         dpc_dbx = P.beta*del_today*P.varphi/(c_tomorrow^P.chic)*...
                        ((dy_dbx(basis_interp,c_tomorrow,pi_tilde_tomorrow)*(pi_tilde_tomorrow-1)*pi_tilde_tomorrow) + y_tomorrow*(2*pi_tilde_tomorrow*dpitilde_dbx(basis_interp) - dpitilde_dbx(basis_interp)));

         % Derivatives with respect to c coefficients
         dpc_dcx = 0*basis_interp';

         % RHS derivative of Value (for Jacobian)
         % Derivatives with respect to a coefficients
         dval_dax  = 0*basis_interp';

         % Derivatives with respect to b coefficients
         dval_dbx  = 0*basis_interp';

         % Derivatives with respect to c coefficients
         dval_dcx  = P.beta*del_today*dv_dcx(basis_interp);

         % Compute GH Integration
         exp_ee_int = exp_ee_int + pi^(-0.5)*G.e_weight(j)*(exp_ee);
         exp_pc_int = exp_pc_int + pi^(-0.5)*G.e_weight(j)*(exp_pc);
         exp_v_int = exp_v_int + pi^(-0.5)*G.e_weight(j)*(exp_v);

         dee_dax_int = dee_dax_int + pi^(-0.5)*G.e_weight(j)*dee_dax;
         dee_dbx_int = dee_dbx_int + pi^(-0.5)*G.e_weight(j)*dee_dbx;
         dee_dcx_int = dee_dcx_int + pi^(-0.5)*G.e_weight(j)*dee_dcx;

         dpc_dax_int = dpc_dax_int + pi^(-0.5)*G.e_weight(j)*dpc_dax;
         dpc_dbx_int = dpc_dbx_int + pi^(-0.5)*G.e_weight(j)*dpc_dbx;
         dpc_dcx_int = dpc_dcx_int + pi^(-0.5)*G.e_weight(j)*dpc_dcx;

         dv_dax_int  = dv_dax_int + pi^(-0.5)*G.e_weight(j)*dval_dax;
         dv_dbx_int  = dv_dbx_int + pi^(-0.5)*G.e_weight(j)*dval_dbx;
         dv_dcx_int  = dv_dcx_int + pi^(-0.5)*G.e_weight(j)*dval_dcx;
     end

     % Residuals
     res1(i) = c_today^(-P.chic) - P.beta*del_today*r_today*exp_ee_int;
     res2(i) = (y_today/c_today^(P.chic))*(P.varphi*(pi_tilde_today-1)*pi_tilde_today-(1-P.theta)-P.theta*(1-P.tau)*w_today) - P.beta*del_today*exp_pc_int;
     if P.chic == 1
        res3(i) = v_today - (log(c_today) - n_today^(P.chin+1)/(P.chin+1)) - P.beta*del_today*exp_v_int;
     else
        res3(i) = v_today - (c_today^(1-P.chic)/(1-P.chic) - n_today^(P.chin+1)/(P.chin+1)) - P.beta*del_today*exp_v_int;
     end

     % Jacobian
     % Jacobian enteries relating to EE
     lhs_jacob_ee(i,:)  = [-P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i)), dc_dbx(C.basis(:,i)),dc_dcx(C.basis(:,i))];
     rhs_jacob_ee(i,:)  = [dee_dax_int, dee_dbx_int,dee_dcx_int];

     % Jacobian enteries relating to PC
     lhs_jacob_pc_ax = (dy_dax(C.basis(:,i),pi_tilde_today)*c_today^(-P.chic) - P.chic*c_today^(-P.chic-1)*dc_dax(C.basis(:,i))*y_today)* ...
                      (P.varphi*(pi_tilde_today - 1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
     lhs_jacob_pc_ax = lhs_jacob_pc_ax + -P.theta*(1-P.tau)*dw_dax(C.basis(:,i),c_today,pi_tilde_today,n_today)* ...
                      y_today/(c_today^(P.chic));

     lhs_jacob_pc_bx = dy_dbx(C.basis(:,i),c_today,pi_tilde_today)/c_today^(P.chic)*(P.varphi*(pi_tilde_today-1)*pi_tilde_today - (1-P.theta) - P.theta*(1-P.tau)*w_today);
     lhs_jacob_pc_bx = lhs_jacob_pc_bx + y_today/(c_today^(P.chic))*(P.varphi*(2*pi_tilde_today*dpitilde_dbx(C.basis(:,i))-dpitilde_dbx(C.basis(:,i))) - P.theta*(1-P.tau)*dw_dbx(C.basis(:,i),c_today,pi_tilde_today,n_today));

     lhs_jacob_pc(i,:)  = [lhs_jacob_pc_ax, lhs_jacob_pc_bx, dpi_dcx(C.basis(:,i))];
     rhs_jacob_pc(i,:)  = [dpc_dax_int,dpc_dbx_int,dpc_dcx_int];

     lhs_jacob_v(i,:) = [dv_dax(C.basis(:,i)), dv_dbx(C.basis(:,i)), dv_dcx(C.basis(:,i))];
     if P.chic == 1
         rhs_jacob_v(i,:) = [c_today^(-1)*dc_dax(C.basis(:,i)) - n_today*dn_dax(C.basis(:,i),pi_tilde_today)+ dv_dax_int, - n_today*dn_dbx(C.basis(:,i),c_today,pi_tilde_today)+dv_dbx_int,dv_dcx_int];
     else
         rhs_jacob_v(i,:) = [c_today*dc_dax(C.basis(:,i)) - n_today*dn_dax(C.basis(:,i),pi_tilde_today)+ dv_dax_int, - n_today*dn_dbx(C.basis(:,i),c_today,pi_tilde_today)+dv_dbx_int,dv_dcx_int];
     end
end

R = [res1,res2,res3];
J = [lhs_jacob_ee-rhs_jacob_ee;...
     lhs_jacob_pc-rhs_jacob_pc;...
     lhs_jacob_v-rhs_jacob_v];

     

     
     
     

