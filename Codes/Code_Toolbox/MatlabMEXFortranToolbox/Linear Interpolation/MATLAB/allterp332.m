function [o1,o2,o3] = allterp332(  x1, x2, x3, ...
                                   x1i, x2i, x3i, ...
                                   pf1 ,pf2, pf3)

% allterp332 linear inter/extrapolation (3 states, 3 policies, 2 stoch comp)
% Inputs:
%   x*      :   Grid
%   x*i     :   Point to evaluate
%   pf*     :   Policy function
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension x4ipts

% Grid lengths
nx1 = length(x1);
nx2 = length(x2);
nx3 = length(x3);

% Number of stochastic realizations
x2ipts = length(x2i);
x3ipts = length(x3i);

% Preallocate output
o1 = zeros(x2ipts,x3ipts);
o2 = zeros(x2ipts,x3ipts);
o3 = zeros(x2ipts,x3ipts);

for i3 = 1:x3ipts
    for i2 = 1:x2ipts 
        s1 = x1(2) - x1(1);
        x1i_min = x1i - x1(1);                     
        loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));

        s2 = x2(2) - x2(1);
        x2i_min = x2i(i2) - x2(1);
        loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));

        s3 = x3(2) - x3(1);
        x3i_min = x3i(i3) - x3(1);
        loc3 = min(nx3-1,max(1,floor(x3i_min/s3) + 1));
        
        xi = [x1i x2i(i2) x3i(i3)];
        xi_left = [x1(loc1) x2(loc2) x3(loc3)];
        xi_right = [x1(loc1+1) x2(loc2+1) x3(loc3+1)];

        w_2 = (xi - xi_left)./ (xi_right - xi_left);
        w_1 = 1 - w_2;
        w1 = [w_1(1) w_2(1)];
        w2 = [w_1(2) w_2(2)];
        w3 = [w_1(3) w_2(3)];

        for m3 = 0:1
            for m2 = 0:1
                for m1 = 0:1
                    o1(i2,i3) = o1(i2,i3) + ...
                            w1(m1+1)*w2(m2+1)*w3(m3+1)*...
                            pf1(loc1+m1,loc2+m2,loc3+m3);
                    o2(i2,i3) = o2(i2,i3) + ...
                            w1(m1+1)*w2(m2+1)*w3(m3+1)*...
                            pf2(loc1+m1,loc2+m2,loc3+m3);

                    o3(i2,i3) = o3(i2,i3) + ...
                            w1(m1+1)*w2(m2+1)*w3(m3+1)*...
                            pf3(loc1+m1,loc2+m2,loc3+m3);
                end
            end
        end        
    end
end