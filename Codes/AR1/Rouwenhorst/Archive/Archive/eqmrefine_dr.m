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

sol = A\b;
i = sol(2*numpts+1:3*numpts);

inx = find(i < 0);
if length(inx) > 1
    inx = inx(end);
end


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
            A1(i,j) = 1 - R.trans_mat(i,j);
            A2(i,j) = -cSIGMA*R.trans_mat(i,j);
            A3(i,j) = cSIGMA;
            A4(i,j) = -cKAPPA;
            if i < inx
                A8(i,j) = -cPHIpi;
            else
                A8(i,j) = 0;
            end
            A9(i,j) = 1;
        else
            A1(i,j) = - R.trans_mat(i,j);
            A2(i,j) = -cSIGMA*R.trans_mat(i,j);
        end
    end
    b1(i) = cRstar - R.del_grid(i);
    b2(i) = 0;
    if i < inx
        b3(i) = cRstar;
    else
        b3(i) = 0;
    end
end

A5 = cBET.*A1;
A7 = A6;

A_up = [A1,A2,A3;...
        A4,A5,A6;...
        A7,A8,A9];

b_up = [b1;b2;b3];

out_up = A_up\b_up;