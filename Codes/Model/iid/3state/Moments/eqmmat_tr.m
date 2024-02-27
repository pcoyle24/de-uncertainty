function [A, b, out] = eqmmat_tr(params,R,numpts)
% [A b] = eqmmat(params,R,numpts)
% Construct Matrix and vector used to solve for solutions of model
%   Inputs:
%       params: a vector of parmaters
%       R: Structure related to Rowenhourst Approx Method
%       numpts: number of elements in discrictized shock vector
%   Outputs:
%       A: Matrix of coefficients
%       b: Vector of coefficients
%       out: Vector of solutions

cBET    = params(1);
cSIGMA  = params(2);
cKAPPA  = params(3);
cPHIpi  = params(4);
cRstar  = params(5);

% Build out in 9 separate blocks of numptsXnumpts.
% Allocate Space for blocks
A1 = zeros(numpts);
A2 = zeros(numpts);
A3 = zeros(numpts);
A4 = zeros(numpts);
A6 = zeros(numpts);
A8 = zeros(numpts);
A9 = zeros(numpts);

b1 = zeros(numpts,1);
b2 = zeros(numpts,1);
b3 = zeros(numpts,1);

for i = 1:numpts
    for j = 1:numpts
        if i == j
            A1(i,j) = 1 - R.s(j);
            A2(i,j) = -cSIGMA*R.s(j);
            A3(i,j) = cSIGMA;
            A4(i,j) = -cKAPPA;
            A8(i,j) = -cPHIpi;
            A9(i,j) = 1;
        else
            A1(i,j) = - R.s(j);
            A2(i,j) = -cSIGMA*R.s(j);
        end
    end
    b1(i) = cRstar - R.del_grid(i);
    b2(i) = 0;
    b3(i) = cRstar;
end

A5 = cBET.*A1;
A7 = A6;

A = [A1,A2,A3;...
    A4, A5,A6;...
    A7,A8,A9];

b = [b1;b2;b3];

out = A\b;


