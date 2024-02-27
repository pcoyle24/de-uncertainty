function [E,sig,Pelb] = sim_pfs_zlb(P,O,C,S)

% E_val = sim_pfs_zlb(P,O,S,G,C,pf)
%   Simulates the model to find Expected Value and Standard Dev of C, Pi, R
% Output:
%   E : Expected Value of C, Pi, R (Structure)
%   sig : Standard Dev of C, Pi, R (Structure)
%   Pelb: Probability of ELB binding

sim = 5001000;
burn = 1000;

%% Calculate Expected Values
disp('Simulating model to calcualte Expected Value')
sum_c_s = 0;
sum_c_d = 0;
sum_pi_s = 0;
sum_pi_d = 0;
sum_r_s = 0;
sum_r_d = 0;
sum_prob_elb_s = 0;
sum_prob_elb_d = 0;

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

r_s = Fallcheb111(O.delbound,sim,del_today,O.n1,C.Ar_s,C.max,C.T,C.P);
r_d = Fallcheb111(O.delbound,sim,del_today,O.n1,C.Ar_d,C.max,C.T,C.P);

parfor i = burn:sim
    if r_s(i) >= 1
        c_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_s,C.max,C.T,C.P);
        pi_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_s,C.max,C.T,C.P);
        r_s_today = r_s(i);
        prob_elb_s =  0;
    else
        c_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_zlb_s,C.max,C.T,C.P);
        pi_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_zlb_s,C.max,C.T,C.P);
        r_s_today = 1;
        prob_elb_s =  1;
    end    
%     sum_c_s = sum_c_s + 100*(c_s_today-S.c_s)/S.c_s;
    sum_c_s = sum_c_s + 100*(c_s_today-1);
    sum_pi_s = sum_pi_s + 400*(pi_s_today-1);
    sum_r_s = sum_r_s + 400*(r_s_today-1);
    sum_prob_elb_s = sum_prob_elb_s + prob_elb_s;
 
    if r_d(i) >= 1
        c_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_d,C.max,C.T,C.P);
        pi_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_d,C.max,C.T,C.P);
        r_d_today = r_d(i);
        prob_elb_d = 0;
    else
        c_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_zlb_d,C.max,C.T,C.P);
        pi_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_zlb_d,C.max,C.T,C.P);
        r_d_today = 1;
        prob_elb_d = 1;
    end    
%     sum_c_d = sum_c_d + 100*(c_d_today-S.c_s)/S.c_s;
    sum_c_d = sum_c_d + 100*(c_d_today-1);
    sum_pi_d = sum_pi_d + 400*(pi_d_today-1);
    sum_r_d = sum_r_d + 400*(r_d_today-1);
    sum_prob_elb_d = sum_prob_elb_d + prob_elb_d;

end


E.c_s = sum_c_s / (sim - burn);
E.pi_s = sum_pi_s / (sim - burn);
E.r_s = sum_r_s / (sim - burn);
Pelb_s = sum_prob_elb_s/ (sim - burn);

E.c_d = sum_c_d / (sim - burn);
E.pi_d = sum_pi_d / (sim - burn);
E.r_d = sum_r_d / (sim - burn);
Pelb_d = sum_prob_elb_d/ (sim - burn);

uncond_prob_s = (1-P.Pd)/(2-P.Ps-P.Pd);
uncond_prob_d = (1-P.Ps)/(2-P.Ps-P.Pd);

E.c = uncond_prob_s*E.c_s + uncond_prob_d*E.c_d;
E.pi = uncond_prob_s*E.pi_s + uncond_prob_d*E.pi_d;
E.r = uncond_prob_s*E.r_s + uncond_prob_d*E.r_d;
Pelb = uncond_prob_s*Pelb_s + uncond_prob_d*Pelb_d;
disp('Expected Value calculated')

%% Calculate Standard Deviation
disp('Simulating model to calcualte Standard Deviation')
sum_c_s_2 = 0;
sum_c_d_2 = 0;
sum_pi_s_2 = 0;
sum_pi_d_2 = 0;
sum_r_s_2 = 0;
sum_r_d_2 = 0;

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

r_s = Fallcheb111(O.delbound,sim,del_today,O.n1,C.Ar_s,C.max,C.T,C.P);
r_d = Fallcheb111(O.delbound,sim,del_today,O.n1,C.Ar_d,C.max,C.T,C.P);

parfor i = burn:sim
    if r_s(i) >= 1
        c_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_s,C.max,C.T,C.P);
        pi_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_s,C.max,C.T,C.P);
        r_s_today = r_s(i);
    else
        c_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_zlb_s,C.max,C.T,C.P);
        pi_s_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_zlb_s,C.max,C.T,C.P);
        r_s_today = 1;
    end    
%     sum_c_s_2 = sum_c_s_2 + (100*(c_s_today-S.c_s)/S.c_s - E.c_s)^2;
    sum_c_s_2 = sum_c_s_2 + (100*(c_s_today-1) - E.c_s)^2;
    sum_pi_s_2 = sum_pi_s_2 + (400*(pi_s_today-1) - E.pi_s)^2;
    sum_r_s_2 = sum_r_s_2 + (400*(r_s_today-1) - E.r_s)^2;
 
    if r_d(i) >= 1
        c_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_d,C.max,C.T,C.P);
        pi_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_d,C.max,C.T,C.P);
        r_d_today = r_d(i);
    else
        c_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ac_zlb_d,C.max,C.T,C.P);
        pi_d_today = Fallcheb111(O.delbound,1,del_today(i),O.n1,C.Ainf_zlb_d,C.max,C.T,C.P);
        r_d_today = 1;
    end    
%     sum_c_d_2 = sum_c_d_2 + (100*(c_d_today-S.c_s)/S.c_s - E.c_d)^2;
    sum_c_d_2 = sum_c_d_2 + (100*(c_d_today-1) - E.c_d)^2;
    sum_pi_d_2 = sum_pi_d_2 + (400*(pi_d_today-1) - E.pi_d)^2;
    sum_r_d_2 = sum_r_d_2 + (400*(r_d_today-1) - E.r_d)^2;
end

sig.c_s = (sum_c_s_2 / (sim - burn))^(1/2);
sig.pi_s = (sum_pi_s_2 / (sim - burn))^(1/2);
sig.r_s = (sum_r_s_2 / (sim - burn))^(1/2);

sig.c_d = (sum_c_d_2 / (sim - burn))^(1/2);
sig.pi_d = (sum_pi_d_2 / (sim - burn))^(1/2);
sig.r_d = (sum_r_d_2 / (sim - burn))^(1/2);

uncond_prob_s = (1-P.Pd)/(2-P.Ps-P.Pd);
uncond_prob_d = (1-P.Ps)/(2-P.Ps-P.Pd);

sig.c = uncond_prob_s*sig.c_s + uncond_prob_d*sig.c_d;
sig.pi = uncond_prob_s*sig.pi_s + uncond_prob_d*sig.pi_d;
sig.r = uncond_prob_s*sig.r_s + uncond_prob_d*sig.r_d;
disp('Standard Deviation calculated')
