function R = eqm_val(weights,P,S,G,O,C,pf)

% [ = eqm_jacob(weights,P,S,G,O,C)
%   Outputs residuals of the equilibrium system of 
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


% Allocate Space for Residuals

res1        = zeros(G.nodes,1);
res1_zlb        = zeros(G.nodes,1);

% Map Coefficients 
C_Av    = weights(1:O.n1+1,1);
C_Av_zlb    = weights(1:O.n1+1,2);

% Get Policy Functions
pf_v =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Av,C.max,C.T,C.P);    
pf_v_zlb =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Av_zlb,C.max,C.T,C.P); 

% Get r_today's cheb weights
C_Ar = Fchebweights11(O.n1,O.del_pts,pf.r,C.T,C.P,C.X);

for i = 1:G.nodes
     del_today         = G.del_gr(i);
     pf_r =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Ar,C.max,C.T,C.P);    
     
     if pf_r(i)  >= 1
         v_today        = pf_v(i);
         
         % Get c_today and n_today
         c_today = pf.c(i);
         n_today = pf.n(i);
              
         % Get Tomorrow's Shock 
         del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

         % Construct GH Integration
        exp_v_int = 0;
        for j = 1:length(G.e_weight)
            r_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar,C.max,C.T,C.P);
            if r_tomorrow >= 1
                v_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av,C.max,C.T,C.P);
                
                % Compute GH Integration
                exp_v_int = exp_v_int + pi^(-0.5)*G.e_weight(j)*v_tomorrow;
            else
                v_tomorrow_zlb = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb,C.max,C.T,C.P);
                
                % Compute GH Integration
                exp_v_int = exp_v_int + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb;
            end
        end
         
        % Residual
        if P.chic == 1
            res1(i) = v_today - log(c_today) + n_today^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_int;
        else
            res1(i) = v_today - c_today^(1-P.chic)/(1-P.chic) + n_today^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_int;
        end                     
     else
        v_today_zlb = pf_v_zlb(i);
        
        % Get c_today and n_today
        c_today_zlb = pf.c_zlb(i);
        n_today_zlb = pf.n_zlb(i);
              
         % Define tomorrow's state variables
        del_tomorrow = P.rho*(del_today-1)+G.e_nodes+1; 
        
        % Construct GH Integration
        exp_v_zlb_int = 0;
        for j = 1:length(G.e_weight)
            r_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar,C.max,C.T,C.P);
            if r_tomorrow >= 1
                v_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av,C.max,C.T,C.P);
                
                % Compute GH Integration
                exp_v_zlb_int = exp_v_zlb_int + pi^(-0.5)*G.e_weight(j)*v_tomorrow;
            else
                v_tomorrow_zlb = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb,C.max,C.T,C.P);
                
                % Compute GH Integration
                exp_v_zlb_int = exp_v_zlb_int + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb;
            end
        end
        
        % Residual
        if P.chic == 1
            res1_zlb(i) = v_today_zlb - log(c_today_zlb) + n_today_zlb^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_zlb_int;
        else
            res1_zlb(i) = v_today_zlb- c_today_zlb^(1-P.chic)/(1-P.chic) + n_today_zlb^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_zlb_int;
        end
     end
 end
 
 R = [res1,res1_zlb];