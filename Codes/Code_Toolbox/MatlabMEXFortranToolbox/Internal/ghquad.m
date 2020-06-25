function [z,w] = ghquad(n)

% ghquad computes the zeros and weights of Gauss-Hermite quadrature formula
% Inputs:
%   n       :   Order of Hermite polynomial
% Outputs:
%   z       :   Zeros of the nth order Hermite polynomial
%   w       :   Weights in Gauss-Hermite quadrature formula

% Coefficients of Hermite polynomial
H = HermitePoly(n+1);
% Corresponding powers
P = (n+2) - (1:n+2);
% Find roots of polynomial
z = sort(roots(HermitePoly(n)));
% Solve for weights as function of roots
w = zeros(n,1);
num = ((2^(n+1))*factorial(n)*(pi^(.5)));
for i = 1:n
    den = sum(H.*(z(i).^P))^2;
    w(i) = num/den;
end