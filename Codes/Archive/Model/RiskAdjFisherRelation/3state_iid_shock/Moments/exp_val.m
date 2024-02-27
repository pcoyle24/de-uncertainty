function [Ev,Eelb] = exp_val(params,pf,R,simval,burn,elb)
rng('default')
cSIGMAd  = params(6);

pf_sum = zeros(size(pf,2),1);
pf_out = zeros(size(pf,2),1);

if nargin > 5;
    elb_sum = 0;
end
u1 = normrnd(0,cSIGMAd/3,simval+burn,1);

for i = 1:simval + burn
    del_today = u1(i);

    if i > burn
        for j = 1:size(pf,2)
            pf_out(j) = allinterp1(R.del_grid',del_today,pf(:,j));
            pf_sum(j) = pf_sum(j) + pf_out(j);
        end
        if nargin > 5
            if pf_out(end) == 0
                elb_sum = elb_sum + 1;
            end
        end
    end
        
end

Ev = pf_sum/simval;

if nargin > 5
    Eelb = elb_sum/simval;
end