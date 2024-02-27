function [o1,o2,o3] = allcheb531_T( z1bnd,z2bnd,z3bnd,z4bnd,z5bnd, ...
                                    z1i,z2i,z3i,z4i,z5i, ...
                                    A1,A2,A3,T,P)

% [o1,o2,o3] = allcheb531_T( z1bnd,z2bnd,z3bnd,z4bnd,z5bnd, ...
%                            z1i,z2i,z3i,z4i,z5i, ...
%                            A1,A2,A3,T,P)
%   Chebyshev inter/extrapolation (5 states, 3 policies, 1 stoch comp)
% Inputs:
%   z*bnd : Min and max of interval on dimension *
%   z*i   : Point to evaluate on dimension *
%   A*    : Coefficients of continuous least squares approximation
%   T     : Coefficients of Chebyshev polynomials
%   P     : Powers of Chebyshev polynomials
% Outputs:
%   o*    : Interpolated/extrapolated values of dimension (x*ipts)

% Determine polynomial degrees
[nn1,nn2,nn3,nn4,nn5] = size(A1);

% Number of stochastic realizations
x5ipts = length(z5i);

% Preallocate output vectors
o1 = zeros(x5ipts,1);
o2 = zeros(x5ipts,1);
o3 = zeros(x5ipts,1);

% Transform point from [a,b] to [-1,1] domain
x1i = 2*(z1i-z1bnd(1))/(z1bnd(2) - z1bnd(1)) - 1;
x2i = 2*(z2i-z2bnd(1))/(z2bnd(2) - z2bnd(1)) - 1;
x3i = 2*(z3i-z3bnd(1))/(z3bnd(2) - z3bnd(1)) - 1;
x4i = 2*(z4i-z4bnd(1))/(z4bnd(2) - z4bnd(1)) - 1;
x5i = 2*(z5i-z5bnd(1))/(z5bnd(2) - z5bnd(1)) - 1;

% Evaluate Chebyshev polynomials
vec1 = sum(T(1:nn1,:).*(x1i).^P(1:nn1,:),2);
vec2 = sum(T(1:nn2,:).*(x2i).^P(1:nn2,:),2);
vec3 = sum(T(1:nn3,:).*(x3i).^P(1:nn3,:),2);
vec4 = sum(T(1:nn4,:).*(x4i).^P(1:nn4,:),2);
for ii = 1:x5ipts
    vec5 = sum(T(1:nn5,:).*(x5i(ii)).^P(1:nn5,:),2);
    [t1 t2 t3 t4 t5] = ndgrid(vec1,vec2,vec3,vec4,vec5);
    
    % Evaluate Chebyshev polynomial at (x1i,x2i)
    temp1 = A1.*t1.*t2.*t3.*t4.*t5;
    temp2 = A2.*t1.*t2.*t3.*t4.*t5;
    temp3 = A3.*t1.*t2.*t3.*t4.*t5;   
    o1(ii) = sum(temp1(:));
    o2(ii) = sum(temp2(:)); 
    o3(ii) = sum(temp3(:));  
end