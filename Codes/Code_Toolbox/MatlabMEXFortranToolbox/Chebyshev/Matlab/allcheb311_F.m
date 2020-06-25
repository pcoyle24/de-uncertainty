function o1 = allcheb311_F(  z1bound,z2bound,z3bound, ...
                             z1i,z2i,z3i, ...
                             A1, ...
                             T,P)
% o1 = allcheb311_F(z1bound,z2bound,z1i,z2i,z3i,A1,T,P)
%   Chebyshev inter/extrapolation (3 states, 1 policy, 1 stoch comp)
% Inputs:
%   z*bound :   Min and max of interval on dimension *
%   z*i     :   Point to evaluate on dimension *
%   A*      :   Coefficients of continuous least squares approximation
%   x*ipts  :   Number of stochastic realizations
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension (x*ipts)
 
% Number of coefficients
[nn1,nn2,nn3] = size(A1);

% Get dimensions
griddim = size(z1i);
nodes = numel(z1i);

% Preallocate output vectors
o1 = zeros(griddim);

% Transform point from [a,b] to [-1,1] domaiparn
x1i = 2*(z1i-z1bound(1))/(z1bound(2) - z1bound(1)) - 1;
x2i = 2*(z2i-z2bound(1))/(z2bound(2) - z2bound(1)) - 1;
x3i = 2*(z3i-z3bound(1))/(z3bound(2) - z3bound(1)) - 1;

% Evaluate Chebyshev polynomials
parfor ii = 1:nodes
    vec1 = sum(T(1:nn1,:).*(x1i(ii)).^P(1:nn1,:),2);
    vec2 = sum(T(1:nn2,:).*(x2i(ii)).^P(1:nn2,:),2);
    vec3 = sum(T(1:nn3,:).*(x3i(ii)).^P(1:nn3,:),2);
    [t1 t2 t3] = ndgrid(vec1,vec2,vec3);
    
    % Evaluate Chebyshev polynomial at (x1i,x2i)
    temp1 = A1.*t1.*t2.*t3;
    o1(ii) = sum(temp1(:));    
end