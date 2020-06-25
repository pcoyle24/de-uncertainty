function [A1 A2 A3] = chebweights53( n1,n2,n3,n4,n5, ...
                                     f1,f2,f3,...
                                     T,P,X)

% [A1 A2 A3] = chebweights53( n1,n2,n3,n4,n5, ...
%                             f1,f2,f3,f4,...
%                             T,P,X)
% Computes the Chebychev polynomial parameters for 3 functions with 5
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
[m1 m2 m3 m4 m5] = size(f1);
% Add 1 to each degree (to account for 0th degree)
nn1 = n1 + 1;
nn2 = n2 + 1;
nn3 = n3 + 1;
nn4 = n4 + 1;
nn5 = n5 + 1;
% Use zeros from desired order of approximation for grid
x1_grid = X(m1,1:m1);
x2_grid = X(m2,1:m2);
x3_grid = X(m3,1:m3);
x4_grid = X(m4,1:m4);
x5_grid = X(m5,1:m5);
[x1_gr x2_gr x3_gr x4_gr x5_gr] = ndgrid(x1_grid,x2_grid,x3_grid,x4_grid,x5_grid);
% Weights to use when calculating coefficients
a1 = zeros(nn1,1);
a1(1) = 1/m1;
a1(2:end) = 2/m1;
a2 = zeros(nn2,1);
a2(1) = 1/m2;
a2(2:end) = 2/m2;
a3 = zeros(nn3,1);
a3(1) = 1/m3;
a3(2:end) = 2/m3;
a4 = zeros(nn4,1);
a4(1) = 1/m4;
a4(2:end) = 2/m4;
a5 = zeros(nn5,1);
a5(1) = 1/m5;
a5(2:end) = 2/m5;
% Calculate coefficients of continuous least squares approx.
A1 = zeros(nn1,nn2,nn3,nn4,nn5);
A2 = A1;
A3 = A1;
for j5 = 1:nn5
    for j4 = 1:nn4
        for j3 = 1:nn3
            for j2 = 1:nn2
                for j1 = 1:nn1
                    nestsum1 = 0;
                    nestsum2 = 0;
                    nestsum3 = 0;
                    T1 = T(j1,:);
                    T2 = T(j2,:);
                    T3 = T(j3,:);
                    T4 = T(j4,:);
                    T5 = T(j5,:);
                    P1 = P(j1,:);
                    P2 = P(j2,:);
                    P3 = P(j3,:);
                    P4 = P(j4,:);
                    P5 = P(j5,:);
                    for kk = 1:nodes
                        % Evaluate Chebyshev polynomials
                        t1 = sum(T1.*(x1_gr(kk)).^P1);
                        t2 = sum(T2.*(x2_gr(kk)).^P2);
                        t3 = sum(T3.*(x3_gr(kk)).^P3);
                        t4 = sum(T4.*(x4_gr(kk)).^P4);
                        t5 = sum(T5.*(x5_gr(kk)).^P5);
                        temp = t1*t2*t3*t4*t5;
                        nestsum1 = nestsum1 + f1(kk)*temp;
                        nestsum2 = nestsum2 + f2(kk)*temp;
                        nestsum3 = nestsum3 + f3(kk)*temp;
                    end
                    temp = a1(j1)*a2(j2)*a3(j3)*a4(j4)*a5(j5);
                    A1(j1,j2,j3,j4,j5) = temp*nestsum1;
                    A2(j1,j2,j3,j4,j5) = temp*nestsum2;
                    A3(j1,j2,j3,j4,j5) = temp*nestsum3;
                end
            end
        end
    end
end
