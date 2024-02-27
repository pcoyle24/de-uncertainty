function [o1] = allinterp1(x1, x1i, pf1) % 
% [o1,o2,o3] = allinterp(x1, x1i, pf1)
% Inputs:
%   x*  : Grids
%   x*i : Point to evaluate
%   pf* : Policy function
% Outputs:
%   o* : Interpolated/extrapolated values of dimension x*ipts

% Grid lengths
nx1 = length(x1);

% Number of stochastic realizations
x1ipts = length(x1i);

% Preallocate output
o1 = zeros(x1ipts,1);
% o2 = zeros(x1ipts,1);
% o3 = zeros(x1ipts,1);

for i1 = 1:x1ipts
    
    s1 = x1(2) - x1(1);
    x1i_min = x1i(i1) - x1(1);
    loc1 = min(nx1 - 1, max(1,floor(x1i_min/s1) + 1));
    
    xi = x1i(i1);
    xi_left = x1(loc1);
    xi_right = x1(loc1+1);
    
    w_2 = (xi - xi_left)./(xi_right - xi_left);
    w_1 = 1 - w_2;
    w1 = [w_1 w_2];
    
    for m1 = 0:1
        o1(i1) = o1(i1) + w1(m1+1)*pf1(loc1+m1);
%         o2(i1) = o2(i1) + w1(m1+1)*pf2(loc1+m1);
%         o3(i1) = o3(i1) + w1(m1+1)*pf3(loc1+m1);
    end
    
end

end

