function R = eqm_val(weights,P,S,G,O,C,pf,CondProb)

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

res1_s        = zeros(G.nodes,1);
res1_zlb_s    = zeros(G.nodes,1);
res1_d        = zeros(G.nodes,1);
res1_zlb_d    = zeros(G.nodes,1);

% Map Coefficients 
C_Av_s      = weights(1:O.n1+1,1);
C_Av_zlb_s  = weights(1:O.n1+1,2);
C_Av_d      = weights(1:O.n1+1,3);
C_Av_zlb_d  = weights(1:O.n1+1,4);

% Get Policy Functions
pf_v_s =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Av_s,C.max,C.T,C.P);    
pf_v_zlb_s =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Av_zlb_s,C.max,C.T,C.P); 
pf_v_d =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Av_d,C.max,C.T,C.P);    
pf_v_zlb_d =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C_Av_zlb_d,C.max,C.T,C.P); 


% Get r_today's cheb weights
C_Ar_s = Fchebweights11(O.n1,O.del_pts,pf.r_s,C.T,C.P,C.X);
C_Ar_d = Fchebweights11(O.n1,O.del_pts,pf.r_d,C.T,C.P,C.X);

for i = 1:G.nodes
     del_today   = G.del_grid(i);
     for k = 1:G.sun_nodes
         cond_prob_sun = CondProb(k,:);
         if k == 1 
             %% Standard Regime
             r_temp =  Fallcheb111(O.delbound,1,del_today,O.n1,C_Ar_s,C.max,C.T,C.P);    
             %% NZLB
             if r_temp >= 1
                 v_today   = pf_v_s(i);
                 % Get c_today and n_today
                 c_today = pf.c_s(i);
                 n_today = pf.n_s(i);

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                 % Construct GH Integration
                 exp_v_s_temp = 0;
                 exp_v_d_temp = 0;
                 for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_s_temp = exp_v_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_s;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_s_temp = exp_v_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_s;
                             end
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_d_temp = exp_v_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_d;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_d_temp = exp_v_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_d;
                             end
                         end
                     end
                     exp_v_s_int = cond_prob_sun*[exp_v_s_temp;exp_v_d_temp];
                 end
                 % Residuals
                if P.chic == 1
                    res1_s(i) = v_today - log(c_today) + n_today^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_s_int;
                else
                    res1_s(i) = v_today - c_today^(1-P.chic)/(1-P.chic) + n_today^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_s_int;
                end   
             else
                 %% ZLB
                 v_today_zlb   = pf_v_zlb_s(i);
                 % Get c_today and n_today
                 c_today_zlb = pf.c_zlb_s(i);
                 n_today_zlb = pf.n_zlb_s(i);

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                 % Construct GH Integration
                 exp_v_zlb_s_temp = 0;
                 exp_v_zlb_d_temp = 0;
                 for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_s_temp = exp_v_zlb_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_s;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_s_temp = exp_v_zlb_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_s;
                             end
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_d_temp = exp_v_zlb_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_d;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_d_temp = exp_v_zlb_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_d;
                             end
                         end
                     end
                     exp_v_zlb_s_int = cond_prob_sun*[exp_v_zlb_s_temp;exp_v_zlb_d_temp];
                 end
                 % Residuals
                if P.chic == 1
                    res1_zlb_s(i) = v_today_zlb - log(c_today_zlb) + n_today_zlb^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_zlb_s_int;
                else
                    res1_zlb_s(i) = v_today_zlb - c_today_zlb^(1-P.chic)/(1-P.chic) + n_today_zlb^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_zlb_s_int;
                end   
             end
         elseif k == 2
             %% Deflationary Regime
             r_temp =  Fallcheb111(O.delbound,1,del_today,O.n1,C_Ar_d,C.max,C.T,C.P);    
             %% NZLB
             if r_temp >= 1
                 v_today   = pf_v_d(i);
                 % Get c_today and n_today
                 c_today = pf.c_d(i);
                 n_today = pf.n_d(i);

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                 % Construct GH Integration
                 exp_v_s_temp = 0;
                 exp_v_d_temp = 0;
                 for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_s_temp = exp_v_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_s;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_s_temp = exp_v_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_s;
                             end
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_d_temp = exp_v_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_d;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_d_temp = exp_v_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_d;
                             end
                         end
                     end
                     exp_v_d_int = cond_prob_sun*[exp_v_s_temp;exp_v_d_temp];
                 end
                 % Residuals
                if P.chic == 1
                    res1_d(i) = v_today - log(c_today) + n_today^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_d_int;
                else
                    res1_d(i) = v_today - c_today^(1-P.chic)/(1-P.chic) + n_today^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_d_int;
                end   
             else
                 %% ZLB
                 v_today_zlb   = pf_v_zlb_d(i);
                 % Get c_today and n_today
                 c_today_zlb = pf.c_zlb_d(i);
                 n_today_zlb = pf.n_zlb_d(i);

                 % Get Tomorrow's Shock 
                 del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

                 % Construct GH Integration
                 exp_v_zlb_s_temp = 0;
                 exp_v_zlb_d_temp = 0;
                 for j = 1:length(G.e_weight)
                     for s = 1:G.sun_nodes
                         if s == 1
                             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_s,C.max,C.T,C.P);
                             if r_s_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_s_temp = exp_v_zlb_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_s;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_s = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_s,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_s_temp = exp_v_zlb_s_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_s;
                             end
                         elseif s == 2
                             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Ar_d,C.max,C.T,C.P);
                             if r_d_tomorrow >= 1
                                 % Get Interpolated Values
                                 v_tomorrow_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_d_temp = exp_v_zlb_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_d;
                             else
                                 % Get Interpolated Values
                                 v_tomorrow_zlb_d = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C_Av_zlb_d,C.max,C.T,C.P);
                                 % Compute GH Integration
                                 exp_v_zlb_d_temp = exp_v_zlb_d_temp + pi^(-0.5)*G.e_weight(j)*v_tomorrow_zlb_d;
                             end
                         end
                     end
                     exp_v_zlb_d_int = cond_prob_sun*[exp_v_zlb_s_temp;exp_v_zlb_d_temp];
                 end
                 % Residuals
                if P.chic == 1
                    res1_zlb_d(i) = v_today_zlb - log(c_today_zlb) + n_today_zlb^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_zlb_d_int;
                else
                    res1_zlb_d(i) = v_today_zlb - c_today_zlb^(1-P.chic)/(1-P.chic) + n_today_zlb^(1+P.chin)/(1+P.chin) - P.beta*del_today*exp_v_zlb_d_int;
                end   
             end
         end
     end
end
                
 
 R = [res1_s,res1_zlb_s,res1_d,res1_zlb_d];