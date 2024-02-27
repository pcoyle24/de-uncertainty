function [A_up, b_up, out_up] = eqmrefine_dr(params,R,numpts,A,b)
% [A b] = eqmmat(params,R,numpts)
% Construct Matrix and vector used to solve for solutions of model
%   Inputs:
%       params: a vector of parmaters
%       R: Structure related to Rowenhourst Approx Method
%       numpts: number of elements in discrictized shock vector
%       A: Matrix of coefficients
%       b: Vector of coefficients
%   Outputs:
%       A_up: Matrix of coefficients (updated)
%       b_up: Vector of coefficients (updated)
%       out_up: Vector of solutions (updated)

cBET    = params(1);
cSIGMA  = params(2);
cKAPPA  = params(3);
cPHIpi  = params(4);
cRstar  = params(5);
cPItarg = params(8);

sol = A\b;
pi = sol(numpts+1:2*numpts);
i = sol(2*numpts+1:3*numpts);
i_implied = cRstar + cPHIpi*pi;

inx_implied = find(i_implied > 0);
inx_i = find(i > 0);
if length(inx_i) > 1
    inx_i = inx_i(end);
end

if isempty(inx_i)
    inx = inx_implied(1);
else
    inx = find(inx_i < inx_implied,1);
end

% Build out in 9 separate blocks of numptsXnumpts.
% Allocate Space for blocks
A11 = zeros(numpts);
A12 = zeros(numpts);
A13 = zeros(numpts);
A21 = zeros(numpts);
A22 = zeros(numpts);
A23 = zeros(numpts);
A32 = zeros(numpts);
A33 = zeros(numpts);

b1 = zeros(numpts,1);
b2 = zeros(numpts,1);
b3 = zeros(numpts,1);

for i = 1:numpts
    for j = 1:numpts
        if i == j
            A11(i,j) = 1 - R.trans_mat(i,j);
            A12(i,j) = -cSIGMA*R.trans_mat(i,j);
            A13(i,j) = cSIGMA;
            A21(i,j) = -cKAPPA;
            A22(i,j) = 1 - cBET*R.trans_mat(i,j);
            if i <= inx
                A32(i,j) = -cPHIpi;
            else
                A32(i,j) = 0;
            end
            A33(i,j) = 1;
        else
            A11(i,j) = - R.trans_mat(i,j);
            A12(i,j) = -cSIGMA*R.trans_mat(i,j);
            A22(i,j) = - cBET*R.trans_mat(i,j);
        end
    end
    b1(i) = cSIGMA*(cRstar - R.del_grid(i));
    b2(i) = cPItarg*(1 - cBET);
    if i <= inx
        b3(i) = cRstar + cPItarg - cPHIpi*cPItarg;
    else
        b3(i) = 0;
    end
end

A31 = A23;

A_up = [A11,A12,A13;...
        A21,A22,A23;...
        A31,A32,A33];

b_up = [b1;b2;b3];

out_up = A_up\b_up;