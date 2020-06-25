function [o1,o2,o3,o4] = allcheb441_F( z1bnd,z2bnd,z3bnd,z4bnd,...
                                       z1i,z2i,z3i,z4i, ...
                                       A1,A2,A3,A4,T,P)
                           
% [o1,o2,o3,o4] = allcheb441_F( z1bnd,z2bnd,z3bnd,z4bnd, ...
%                               z1i,z2i,z3i,z4i, ...
%                               A1,A2,A3,A4,T,P)
%   Chebyshev inter/extrapolation (4 states, 4 policies, 1 stoch comp)
% Inputs:
%   z*bnd : Min and max of interval on dimension *
%   z*i   : Point to evaluate on dimension *
%   A*    : Coefficients of continuous least squares approximation
%   T     : Coefficients of Chebyshev polynomials
%   P     : Powers of Chebyshev polynomials
% Outputs:
%   o*    : Interpolated/extrapolated values of dimension (x*ipts)

% Determine polynomial degrees
[nn1,nn2,nn3,nn4] = size(A1);

% Get dimensions
griddim = size(z1i);
nodes = numel(z1i);

% Preallocate output vectors
o1 = zeros(griddim);
o2 = o1;
o3 = o1;
o4 = o1;

% Transform point from [a,b] to [-1,1] domain
x1i = 2*(z1i-z1bnd(1))/(z1bnd(2) - z1bnd(1)) - 1;
x2i = 2*(z2i-z2bnd(1))/(z2bnd(2) - z2bnd(1)) - 1;
x3i = 2*(z3i-z3bnd(1))/(z3bnd(2) - z3bnd(1)) - 1;
x4i = 2*(z4i-z4bnd(1))/(z4bnd(2) - z4bnd(1)) - 1;

% Evaluate Chebyshev polynomials
parfor ii = 1:nodes
    vec1 = sum(T(1:nn1,:).*(x1i(ii)).^P(1:nn1,:),2);
    vec2 = sum(T(1:nn2,:).*(x2i(ii)).^P(1:nn2,:),2);
    vec3 = sum(T(1:nn3,:).*(x3i(ii)).^P(1:nn3,:),2);
    vec4 = sum(T(1:nn4,:).*(x4i(ii)).^P(1:nn4,:),2);
    [t1 t2 t3 t4] = ndgrid(vec1,vec2,vec3,vec4);
    
    % Evaluate Chebyshev polynomial at (x1i,x2i)
    temp1 = A1.*t1.*t2.*t3.*t4;
    temp2 = A2.*t1.*t2.*t3.*t4; 
    temp3 = A3.*t1.*t2.*t3.*t4; 
    temp4 = A4.*t1.*t2.*t3.*t4;   
    o1(ii) = sum(temp1(:));
    o2(ii) = sum(temp2(:));  
    o3(ii) = sum(temp3(:));  
    o4(ii) = sum(temp4(:));   
end