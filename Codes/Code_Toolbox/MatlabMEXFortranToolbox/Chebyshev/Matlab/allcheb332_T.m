function [o1 o2 o3] = allcheb332_T( z1bnd,z2bnd,z3bnd, ...
                                    z1i,z2i,z3i, ...
                                    A1,A2,A3,T,P)

% [o1 o2 o3] = allcheb332_T( z1bnd,z2bnd,z3bnd, ...
%                            z1i,z2i,z3i, ...
%                            A1,A2,A3,T,P)
% 	Chebyshev inter/extrapolation (3 states, 1 policies, 1 stoch comp)
% Inputs:
%   z*bound :   Min and max of interval on dimension *
%   z*i     :   Point to evaluate on dimension *
%   A*      :   Coefficients of continuous least squares approximation
%   T       :   Coefficients of Chebyshev polynomials
%   P       :   Powers of Chebyshev polynomials
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension (x*ipts)

% Number of coefficients
[nn1,nn2,nn3] = size(A1);

% Number of stochastic realizations
x2ipts = length(z2i);
x3ipts = length(z3i);

% Preallocate output vectors
o1 = zeros(x2ipts,x3ipts);
o2 = o1;
o3 = o1;

% Transform point from [a,b] to [-1,1] domain
x1i = 2*(z1i-z1bnd(1))/(z1bnd(2) - z1bnd(1)) - 1;
x2i = 2*(z2i-z2bnd(1))/(z2bnd(2) - z2bnd(1)) - 1;
x3i = 2*(z3i-z3bnd(1))/(z3bnd(2) - z3bnd(1)) - 1;

% Evaluate Chebyshev polynomials
vec1 = sum(T(1:nn1,:).*(x1i).^P(1:nn1,:),2);
for i3 = 1:x3ipts
    for i2 = 1:x2ipts
        vec2 = sum(T(1:nn2,:).*(x2i(i2)).^P(1:nn2,:),2);
        vec3 = sum(T(1:nn3,:).*(x3i(i3)).^P(1:nn3,:),2);
        [t1 t2 t3] = ndgrid(vec1,vec2,vec3);
        
        % Evaluate Chebyshev polynomial at (x1i,x2i)
        temp1 = A1.*t1.*t2.*t3; 
        temp2 = A2.*t1.*t2.*t3; 
        temp3 = A3.*t1.*t2.*t3; 
        o1(i2,i3) = sum(temp1(:)); 
        o2(i2,i3) = sum(temp2(:)); 
        o3(i2,i3) = sum(temp3(:));  
    end
end