function [A_up, b_up, out_up] = eqmrefine(params,R,numpts,A,b,pim)
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
cPItarg = params(6);
pi_m = pim;

sol = A\b;
i = sol(2*numpts+1:3*numpts);

inx = find(i < 0);
if length(inx) > 1
    inx = inx(end);
end

% Find Middle State
mid = (numpts+1)/2;

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

% Note that we are fixing pi_m. This part of the matrix will account for
% that.
Apim1 = zeros(1,numpts); 
Apim2 = zeros(1,numpts); 
Apim3 = zeros(1,numpts); 

by = zeros(numpts,1);
bpi = zeros(numpts,1);
bi = zeros(numpts,1);

for i = 1:numpts
    for j = 1:numpts
        if i == j
            A11(i,j) = 1 - R.trans_mat(i,j);
            A12(i,j) = -cSIGMA*R.trans_mat(i,j);
            A13(i,j) = cSIGMA;
            A21(i,j) = -cKAPPA;
            A22(i,j) = 1 - cBET*R.trans_mat(i,j);
            if i >= inx
                A32(i,j) = 0;
            else
                A32(i,j) = -cPHIpi;
            end
            A33(i,j) = 1;
            if i == mid
                Apim2(j) = 1;
            end
        else          
            A11(i,j) = -R.trans_mat(i,j);
            A12(i,j) = -cSIGMA*R.trans_mat(i,j);
            A22(i,j) = - cBET*R.trans_mat(i,j);
        end
    end
    by(i) = cSIGMA*(cRstar - R.del_grid(i));
    bpi(i) = cPItarg*(1 - cBET);
    if i >= inx
        bi(i) = 0;
    else
        bi(i) = cRstar + cPItarg - cPHIpi*cPItarg;
    end
end

A31 = A23;

% We consider the 'demeaned Euler Equation' which is y_i - y_mid
% for i ~= mid.
Ay = [A11, A12, A13];
Aymid = repmat(Ay(mid,:),numpts,1);


Ay = Ay - Aymid;
Api = [A21, A22, A23];
Ai = [A31,A32,A33];
Apim = [Apim1,Apim2,Apim3];

bymid = repmat(by(mid),numpts,1);
by = by - bymid;

A_up = [Ay;...
    Api;...
    Ai;
    Apim];

b_up = [by;bpi;bi;pi_m];

% Delete row of zeros.
A_up(mid,:) = [];
b_up(mid,:) = [];

out_up = A_up\b_up;

% Note this is a vecor of 'demeaned EE': y_i - y_mid
y = out_up(1:numpts);
pi = out_up(numpts+1:2*numpts);
i = out_up(2*numpts+1:3*numpts);

Epi = pi'*R.trans_mat(mid,:)';

% Add y_mid back to get y_i
y = y + (pi(mid) - cPItarg - cBET*(Epi - cPItarg))/cKAPPA;

out_up = [y;pi;i];
