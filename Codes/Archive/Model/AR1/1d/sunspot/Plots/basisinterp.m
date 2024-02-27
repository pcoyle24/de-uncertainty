function b = basisinterp(  z1bound, ...
                             z1i, ...
                             A1, ...
                             T,P)
                         
% b = basisinterp(  z1bound,z1i,A1,T,P)
%   Chebyshev Basis inter/extrapolation (1 states, 1 policy, 1 stoch comp)
% Inputs:
%   z*bound :   Min and max of interval on dimension *
%   z*i     :   Point to evaluate on dimension *
%   A*      :   Coefficients of continuous least squares approximation
%   x*ipts  :   Number of stochastic realizations
% Outputs:
%   b*      :   Interpolated/extrapolated values of dimension (x*ipts)                         
                         
% Number of coefficients
[nn1,~] = size(A1);

nodes = numel(z1i);

% Transform point from [a,b] to [-1,1] domain
x1i = 2*(z1i-z1bound(1))/(z1bound(2) - z1bound(1)) - 1;

% Preallocate Basis 
b = zeros(nn1,nodes);

% Evaluate Chebyshev polynomials
for ii = 1:nodes
    vec1 = sum(T(1:nn1,1:nn1).*(x1i(ii)).^P(1:nn1,1:nn1),2);
    t1 = ndgrid(vec1);
    b(:,ii) = t1;
end