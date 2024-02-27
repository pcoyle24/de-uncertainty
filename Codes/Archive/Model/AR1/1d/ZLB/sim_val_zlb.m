function E_val = sim_val_zlb(P,O,C)

% E_val = sim_val(P,O,S,G,C,pf)
%   Simulates the model to find Expected Value
% Output:
%   E_val : Expected Value of Economy

sim = 5001000;
burn = 1000;
sum_v = 0;

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

r_today = Fallcheb111(O.delbound,sim,del_today,O.n1,C.Ar,C.max,C.T,C.P);

parfor i = burn:sim
    if r_today(i) >= 1
        v_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Av,C.max,C.T,C.P);
    else
        v_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Av_zlb,C.max,C.T,C.P);
    end
    
    sum_v = sum_v + v_today;
end

% E_val = mean(v_today);
E_val = sum_v / (sim - burn);

