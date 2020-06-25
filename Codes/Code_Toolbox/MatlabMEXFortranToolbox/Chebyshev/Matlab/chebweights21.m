function A1 = chebweights21( n1,n2,...
                             f1,...
                             T,P,X)

% A1 = chebweights21( n1,n2,...
%                     f1,...
%                     T,P,X)
% Computes the Chebychev polynomial parameters for 1 functions with 2
%   dimensions.
% Inputs:
%     n*    : degree of Chebychev polynomial, must be less than or equal to
%               the corresponding strides of f*
%     f*    : Function to parameterize, all of same dimension
%     T     :   Chebyshev Polynomial Coefficients
%     P     :   Chebyshev Polynomial Powers
%     X     :   Chebyshev Polynomial Zeros
% Outputs:
%     A*    : Least squares coefficients for function *

% Dimensions
nodes = numel(f1);
[m1 m2] = size(f1);
% Add 1 to each degree (to account for 0th degree)
nn1 = n1 + 1;
nn2 = n2 + 1;
% Use zeros from desired order of approximation for grid
x1_grid = X(m1,1:m1);
x2_grid = X(m2,1:m2);
[x1_gr x2_gr] = ndgrid(x1_grid,x2_grid);
% Weights when calculating coefficients
a1 = zeros(nn1,1);
a1(1) = 1/m1;
a1(2:end) = 2/m1;
a2 = zeros(nn2,1);
a2(1) = 1/m2;
a2(2:end) = 2/m2;
% Calculate coefficients of continuous least squares approx.
A1 = zeros(nn1,nn2);
for j2 = 1:nn2
    for j1 = 1:nn1
        nestsum1 = 0;
        T1 = T(j1,:);
        T2 = T(j2,:);
        P1 = P(j1,:);
        P2 = P(j2,:);
        for kk = 1:nodes
            % Evaluate Chebyshev polynomials
            t1 = sum(T1.*(x1_gr(kk)).^P1);
            t2 = sum(T2.*(x2_gr(kk)).^P2);
            nestsum1 = nestsum1 + f1(kk)*t1*t2;
        end
        temp = a1(j1)*a2(j2);
        A1(j1,j2) = temp*nestsum1;
    end
end