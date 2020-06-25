function [o1,o2,o3] = allterp531( x1, x2, x3, x4, x5, ...
                                  x1i, x2i, x3i, x4i, x5i, ...
                                  pf1, pf2, pf3)

% allterp531 linear inter/extrapolation (5 states, 3 policies, 1 stoch comp)
% Inputs:
%   x*      :   Grid
%   x*i     :   Point to evaluate
%   pf*     :   Policy function
%   x*ipts  :   Number of stochastic realizations
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension x*ipts

% Grid lengths
nx1 = length(x1);
nx2 = length(x2);
nx3 = length(x3);
nx4 = length(x4);
nx5 = length(x5);

% Number of stochastic realizations
x5ipts = length(x5i);

% Preallocate output
o1 = zeros(x5ipts,1);
o2 = zeros(x5ipts,1);
o3 = zeros(x5ipts,1);

for i5 = 1:x5ipts
    s1 = x1(2) - x1(1);
    x1i_min = x1i - x1(1);
    loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));

    s2 = x2(2) - x2(1);
    x2i_min = x2i - x2(1);
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));

    s3 = x3(2) - x3(1);
    x3i_min = x3i - x3(1);
    loc3 = min(nx3-1,max(1,floor(x3i_min/s3) + 1));

    s4 = x4(2) - x4(1);
    x4i_min = x4i(i5) - x4(1);
    loc4 = min(nx4-1,max(1,floor(x4i_min/s4) + 1));
    
    s5 = x5(2) - x5(1);
    x5i_min = x5i(i5) - x5(1);
    loc5 = min(nx5-1,max(1,floor(x5i_min/s5) + 1));    
    
    xi = [x1i x2i x3i x4i x5i(i5)];
    xi_left = [x1(loc1) x2(loc2) x3(loc3) x4(loc4) x5(loc5)];
    xi_right = [x1(loc1+1) x2(loc2+1) x3(loc3+1) x4(loc4+1) x5(loc5+1)];

    w_2 = (xi - xi_left)./ (xi_right - xi_left);
    w_1 = 1 - w_2;
    w1 = [w_1(1) w_2(1)];
    w2 = [w_1(2) w_2(2)];
    w3 = [w_1(3) w_2(3)];
    w4 = [w_1(4) w_2(4)];
    w5 = [w_1(5) w_2(5)];
    
    for m5 = 0:1
        for m4 = 0:1
            for m3 = 0:1
                for m2 = 0:1
                    for m1 = 0:1
                        o1(i5) = o1(i5) + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*pf1(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5);
                        o2(i5) = o2(i5) + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*pf2(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5);
                        o3(i5) = o3(i5) + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*pf3(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5);
                    end
                end
            end
        end
    end
end