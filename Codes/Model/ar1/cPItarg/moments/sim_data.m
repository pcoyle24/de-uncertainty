function [pf_sim,elb_sim] = sim_data(params,pf,R,simval,burn,elb)

rng('default')

cRHO    = params(6);
cSIGMAd  = params(7);

del_yesterday = 0;

pf_sim = zeros(simval, 1);

if nargin > 5;
    elb_sim = zeros(simval, 1);
end

for i = 1:simval + burn
    
    u1 = rand;
    u2 = rand;

    ed_today = cSIGMAd*(-2*log(u1))^(1/2)*cos(2*pi*u2);
    del_today = cRHO*(del_yesterday) + ed_today;
    
    if i > burn
        pf_out = allinterp1(R.del_grid',del_today,pf);
        pf_sim(i-burn) =  pf_out;
        if nargin > 5
            if pf_out < 1e-10
                elb_sim(i - burn) = 1;
            end
        end
    end
    
    del_yesterday = del_today;  
end
