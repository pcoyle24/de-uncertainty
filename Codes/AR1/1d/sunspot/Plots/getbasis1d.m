function o1 = getbasis1d(z1bound,z1i,n1,T,P)

% o1 = getbasis(z1bound,z1i,z2i,T,P)
%   Compute Chebyshev basis  (1 states, 1 stoch comp)
% Inputs:
%   z*bound :   Min and max of interval on dimension *
%   z*i     :   Point to evaluate on dimension *
%   nn1     :   Function degree of approximation
%   x*ipts  :   Number of stochastic realizations
% Outputs:
%   o*      :   Interpolated/extrapolated values of dimension (x*ipts)


% Number of coefficients
nn1 = n1+1;

% Get dimensions
griddim = size(z1i);
nodes = numel(z1i);

% Transform point from [a,b] to [-1,1] domain
x1i = 2*(z1i-z1bound(1))/(z1bound(2) - z1bound(1)) - 1;

% Preallocate Basis 
o1 = zeros(nn1,nodes);

% Evaluate Chebyshev polynomials
for ii = 1:nodes
    vec1 = sum(T(1:nn1,:).*(x1i(ii)).^P(1:nn1,:),2);
    t1 = ndgrid(vec1);
    o1(:,ii) = t1;
end