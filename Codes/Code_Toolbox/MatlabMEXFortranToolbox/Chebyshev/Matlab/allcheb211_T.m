function o1 = allcheb211_T( z1bnd,z2bnd, ...
                            z1i,z2i, ...
                            A1,T,P)

% o1 = allcheb211_T(z1bnd,z2bnd, ...
%                   z1i,z2i, ...
%                   A1,T,P)
% 	Chebyshev inter/extrapolation (2 states, 1 policies, 1 stoch comp)
% Inputs:
%   z*bound :   Min and max of interval on dimension *
%   z*i     :   Point to evaluate on dimension *
%   A*      :   Coefficients of continuous least squares approximation
%   T       :   Coefficients of Chebyshev polynomials
%   P       :   Powers of Chebyshev polynomials
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension (x*ipts)

% Number of coefficients
[nn1,nn2] = size(A1);

% Number of stochastic realizations
x2ipts = length(z2i);

% Preallocate output vectors
o1 = zeros(x2ipts,1);

% Transform point from [a,b] to [-1,1] domain
x1i = 2*(z1i-z1bnd(1))/(z1bnd(2) - z1bnd(1)) - 1;
x2i = 2*(z2i-z2bnd(1))/(z2bnd(2) - z2bnd(1)) - 1;

% Evaluate Chebyshev polynomials
vec1 = sum(T(1:nn1,:).*(x1i).^P(1:nn1,:),2);
for ii = 1:x2ipts
    vec2 = sum(T(1:nn2,:).*(x2i(ii)).^P(1:nn2,:),2);
    t1 = vec1(:,ones(nn2,1));
    t2 = vec2(:,ones(nn1,1))';
    
    % Evaluate Chebyshev polynomial at (x1i,x2i)
    temp1 = A1.*t1.*t2; 
    o1(ii) = sum(temp1(:));  
end