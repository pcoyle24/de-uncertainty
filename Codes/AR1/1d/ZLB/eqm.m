function R = eqm(weights,P,S,G,O,C)

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

% Allocate Space for Residuals

res1        = zeros(G.nodes,1);
res2        = zeros(G.nodes,1);
res1_zlb        = zeros(G.nodes,1);
res2_zlb        = zeros(G.nodes,1);

% Map Coefficients 
C.Ac    = weights(1:O.n1+1,1);
C.Ainf  = weights(1:O.n1+1,2);
C.Ac_zlb    = weights(1:O.n1+1,3);
C.Ainf_zlb  = weights(1:O.n1+1,4);

% Get Policy Functions

 pf_c           = allcheb111(O.delbound,G.del_grid,C.Ac,C.T,C.P);    
 pf_inf         = allcheb111(O.delbound,G.del_grid,C.Ainf,C.T,C.P); 
 
 pf_r           = pf_inf.^(P.phi_pi)/P.beta;
 
 pf_c_zlb       = allcheb111(O.delbound,G.del_grid,C.Ac_zlb,C.T,C.P);    
 pf_inf_zlb     = allcheb111(O.delbound,G.del_grid,C.Ainf_zlb,C.T,C.P); 
 
 pf_r_zlb       = ones(length(pf_r),1);    
 
 for i = 1:G.nodes
     del_today         = G.del_gr(i);
     
     c_temp            = pf_c(i);
     inf_temp          = pf_inf(i);
     r_temp            = inf_temp^(P.phi_pi)/P.beta;
     
     if r_temp  >= 1
         c_today       = c_temp;
         inf_today     = inf_temp;
         r_today       = r_temp;
         % Build out other PFs  
         n_today       = (c_today/(1-(P.varphi/2)*(inf_today-1)^2));
         
         % Get Tomorrow's Shock 
         del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;
         
         % Get r_today's cheb weights
         
         C.Ar = chebweights11(O.n1,pf_r,C.T,C.P,C.X);

         %For GH Integration
         exp_ee_int = 0;
         exp_pc_int = 0;

         for j = 1:length(G.e_weight)
             r_tomorrow = allcheb111(O.delbound,del_tomorrow(j),C.Ar,C.T,C.P);
             
             if r_tomorrow >= 1
                 % Get Interpolated Values
                 c_tomorrow         = allcheb111(O.delbound,del_tomorrow(j),C.Ac,C.T,C.P);
                 inf_tomorrow       = allcheb111(O.delbound,del_tomorrow(j),C.Ainf,C.T,C.P);

                 % Build out other PF Values
                 n_tomorrow         = (c_tomorrow/(1-(P.varphi/2)*(inf_tomorrow-1)^2)); 

                 % Compute GH Integration

                 exp_ee             = c_tomorrow^(-P.chic)*inf_tomorrow^(-1);
                 exp_pc             = (n_tomorrow/c_tomorrow^(P.chic))*P.varphi*(inf_tomorrow - 1)*inf_tomorrow;
                 
                 exp_ee_int         = exp_ee_int + pi^(-0.5)*G.e_weight(j)*exp_ee;
                 exp_pc_int         = exp_pc_int + pi^(-0.5)*G.e_weight(j)*exp_pc;
             else 
                 % Get Interpolated Values
                 c_tomorrow_zlb     = allcheb111(O.delbound,del_tomorrow(j),C.Ac_zlb ,C.T,C.P);
                 inf_tomorrow_zlb   = allcheb111(O.delbound,del_tomorrow(j),C.Ainf_zlb ,C.T,C.P);

                 % Build out other PF Values
                 n_tomorrow_zlb     = (c_tomorrow_zlb /(1-(P.varphi/2)*(inf_tomorrow_zlb -1)^2)); 

                 % Compute GH Integration

                 exp_ee_zlb         = c_tomorrow_zlb ^(-P.chic)*inf_tomorrow_zlb ^(-1);
                 exp_pc_zlb         = (n_tomorrow_zlb /c_tomorrow_zlb ^(P.chic))*P.varphi*(inf_tomorrow_zlb  - 1)*inf_tomorrow_zlb;
                 
                 exp_ee_int         = exp_ee_int + pi^(-0.5)*G.e_weight(j)*exp_ee_zlb;
                 exp_pc_int         = exp_pc_int + pi^(-0.5)*G.e_weight(j)*exp_pc_zlb;
             end            
         end
         % Residuals
         res1(i) = c_today^(-P.chic) - P.beta*del_today*r_today*exp_ee_int;
         res2(i) = (n_today/c_today^(P.chic))*(P.varphi*(inf_today-1)*inf_today-(1-P.theta)-P.theta*n_today^(P.chin)*c_today^(P.chic)) - P.beta*del_today*exp_pc_int;
                 
     else
         c_today_zlb        = pf_c_zlb(i);
         inf_today_zlb      = pf_inf_zlb(i);
         r_today_zlb        = S.r_zlb;
         % Build out other PFs  
         n_today_zlb        = (c_today_zlb/(1-(P.varphi/2)*(inf_today_zlb-1)^2));
         
         % Get Tomorrow's Shock 
         del_tomorrow       = P.rho*(del_today - 1) + 1 + G.e_nodes;
         
         % Get r_today's cheb weights
         
         C.Ar_zlb = chebweights11(O.n1,pf_r_zlb,C.T,C.P,C.X);

         %For GH Integration
         exp_ee_int_zlb = 0;
         exp_pc_int_zlb = 0;

         for j = 1:length(G.e_weight)
             r_tomorrow_zlb = allcheb111(O.delbound,del_tomorrow(j),C.Ar,C.T,C.P);
             
             if r_tomorrow_zlb >= 1
                 % Get Interpolated Values
                 c_tomorrow         = allcheb111(O.delbound,del_tomorrow(j),C.Ac,C.T,C.P);
                 inf_tomorrow       = allcheb111(O.delbound,del_tomorrow(j),C.Ainf,C.T,C.P);

                 % Build out other PF Values
                 n_tomorrow         = (c_tomorrow/(1-(P.varphi/2)*(inf_tomorrow-1)^2)); 

                 % Compute GH Integration

                 exp_ee             = c_tomorrow^(-P.chic)*inf_tomorrow^(-1);
                 exp_pc             = (n_tomorrow/c_tomorrow^(P.chic))*P.varphi*(inf_tomorrow - 1)*inf_tomorrow;
                 
                 exp_ee_int_zlb     = exp_ee_int_zlb + pi^(-0.5)*G.e_weight(j)*exp_ee;
                 exp_pc_int_zlb     = exp_pc_int_zlb + pi^(-0.5)*G.e_weight(j)*exp_pc;
             else 
                 % Get Interpolated Values
                 c_tomorrow_zlb     = allcheb111(O.delbound,del_tomorrow(j),C.Ac_zlb ,C.T,C.P);
                 inf_tomorrow_zlb   = allcheb111(O.delbound,del_tomorrow(j),C.Ainf_zlb ,C.T,C.P);

                 % Build out other PF Values
                 n_tomorrow_zlb     = (c_tomorrow_zlb /(1-(P.varphi/2)*(inf_tomorrow_zlb -1)^2)); 

                 % Compute GH Integration

                 exp_ee_zlb         = c_tomorrow_zlb ^(-P.chic)*inf_tomorrow_zlb ^(-1);
                 exp_pc_zlb         = (n_tomorrow_zlb /c_tomorrow_zlb ^(P.chic))*P.varphi*(inf_tomorrow_zlb  - 1)*inf_tomorrow_zlb;
                 
                 exp_ee_int_zlb     = exp_ee_int_zlb + pi^(-0.5)*G.e_weight(j)*exp_ee_zlb;
                 exp_pc_int_zlb     = exp_pc_int_zlb + pi^(-0.5)*G.e_weight(j)*exp_pc_zlb;
             end            
         end
      
         res1_zlb(i) = c_today_zlb^(-P.chic) - P.beta*del_today*r_today_zlb*exp_ee_int_zlb;
         res2_zlb(i) = (n_today_zlb/c_today_zlb^(P.chic))*(P.varphi*(inf_today_zlb-1)*inf_today_zlb-(1-P.theta)-P.theta*n_today_zlb^(P.chin)*c_today_zlb^(P.chic)) - P.beta*del_today*exp_pc_int_zlb;
     end
 end
 
 R = [res1,res2,res1_zlb,res2_zlb];
% disp(num2str(sum(sum(R.^2))));
     
     
     

