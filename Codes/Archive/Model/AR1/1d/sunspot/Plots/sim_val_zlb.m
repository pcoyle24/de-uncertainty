function E_val = sim_val_zlb(P,O,C)

% E_val = sim_val(P,O,S,G,C,pf)
%   Simulates the model to find Expected Value
% Output:
%   E_val : Expected Value of Economy

sim = 5001000;
burn = 1000;
sum_v_s = 0;
sum_v_d = 0;

del_yesterday = 1;
del_today = zeros(sim,1);
rn = rand(sim,2);

for i = 1:sim
    u1 = rn(i,1);
    u2 = rn(i,2);
    
    ed_today = P.sigma*(-2*log(u1))^(1/2)*cos(2*pi*u2);
    del_today(i) = P.rho*(del_yesterday - 1) + 1 + ed_today;    
    del_yesterday = del_today(i); 
end

r_s_today = Fallcheb111(O.delbound,sim,del_today,O.n1,C.Ar_s,C.max,C.T,C.P);
r_d_today = Fallcheb111(O.delbound,sim,del_today,O.n1,C.Ar_d,C.max,C.T,C.P);

parfor i = burn:sim
    if r_s_today(i) >= 1
        v_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Av_s,C.max,C.T,C.P);
    else
        v_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Av_zlb_s,C.max,C.T,C.P);
    end    
    sum_v_s = sum_v_s + v_s_today;
% end
% parfor i = burn:sim
 
    if r_d_today(i) >= 1
        v_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Av_d,C.max,C.T,C.P);
    else
        v_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Av_zlb_d,C.max,C.T,C.P);
    end    
    sum_v_d = sum_v_d + v_d_today;
end

% E_val = mean(v_today);
E_val_s = sum_v_s / (sim - burn);
E_val_d = sum_v_d / (sim - burn);

uncond_prob_s = (1-P.Pd)/(2-P.Ps-P.Pd);
uncond_prob_d = (1-P.Ps)/(2-P.Ps-P.Pd);

E_val = uncond_prob_s*E_val_s + uncond_prob_d*E_val_d;
