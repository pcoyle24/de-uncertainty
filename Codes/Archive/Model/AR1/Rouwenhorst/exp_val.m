function [Ev,Eelb] = exp_val(params,pf,R,simval,burn,elb)

cRHO    = params(6);
cSIGMAd  = params(7);

del_yesterday = 0;

pf_sum = 0;

if nargin > 5;
    elb_sum = 0;
end

for i = 1:simval + burn
    
    u1 = rand;
    u2 = rand;

    ed_today = cSIGMAd*(-2*log(u1))^(1/2)*cos(2*pi*u2);
    del_today = cRHO*(del_yesterday) + ed_today;
    
    if i > burn
        pf_out = allinterp1(R.del_grid',del_today,pf);
        pf_sum = pf_sum + pf_out;
        if nargin > 5
            if pf_out == 0
                elb_sum = elb_sum + 1;
            end
        end
    end
    
    del_yesterday = del_today; 
    
end

Ev = pf_sum/simval;

if nargin > 5
    Eelb = elb_sum/simval;
end