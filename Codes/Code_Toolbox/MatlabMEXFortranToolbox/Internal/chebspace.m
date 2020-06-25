function z_grid = chebspace(a,b,n)

% z_grid = chebspace(a,b,n)
%   Outputs grids with points that correspond to the zeros in the
%   (n-1)th order Chebyshev polynomial.
% Inputs:
%     a     :   Minimum of original interval
%     b     :   Maximum of original interval
%     n     :   Number of grid points
% Output:
%     z_grid:   Vector of zeros in [a,b] domain

% Chebyshev polynomial coefficient, power matrix and zeros
x_grid = zeros(1,n); 
for j = 1:n 
    % Zeros of polynomial
    x_grid(j) = -cos((2*j-1)*pi/(2*n));
end

% Transform x_grid in [-1,1] to z in [a,b] (rescaling)
z_grid = (x_grid+1)*(b - a)/2 + a;
