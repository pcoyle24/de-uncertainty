function [A, b, out] = eqmmat(params,R,numpts,pim)
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
pi_m = pim;

% Find Middle State
mid = (numpts+1)/2;

% Build out in 9 separate blocks of numptsXnumpts.
% Allocate Space for blocks
A1 = zeros(numpts);
A2 = zeros(numpts);
A3 = zeros(numpts);
A4 = zeros(numpts);
A5 = zeros(numpts);
A6 = zeros(numpts);
A8 = zeros(numpts);
A9 = zeros(numpts);

Apim1 = zeros(1,numpts); 
Apim2 = zeros(1,numpts); 
Apim3 = zeros(1,numpts); 

by = zeros(numpts,1);
bpi = zeros(numpts,1);
bi = zeros(numpts,1);

for i = 1:numpts
    for j = 1:numpts
        if i == j
            A1(i,j) = 1 - R.trans_mat(i,j);           
            A2(i,j) = -cSIGMA*R.trans_mat(i,j);
            A3(i,j) = cSIGMA;
            A4(i,j) = -cKAPPA;
            A5(i,j) = cBET*(1 - R.trans_mat(i,j));
            A8(i,j) = -cPHIpi;
            A9(i,j) = 1;
            if i == mid
                Apim2(j) = 1;
            end
        else
            A1(i,j) = -R.trans_mat(i,j);
            A2(i,j) = -cSIGMA*R.trans_mat(i,j);
            A5(i,j) = - cBET*R.trans_mat(i,j);
        end
    end
    by(i) = cRstar - R.del_grid(i);
    bpi(i) = 0;
    bi(i) = cRstar;
end

A7 = A6;

% We consider the 'demeaned Euler Equation' which is y_i - y_mid
% for i ~= mid.
Ay = [A1, A2, A3];
Aymid = repmat(Ay(mid,:),numpts,1);

Ay = Ay - Aymid;
Api = [A4, A5, A6];
Ai = [A7,A8,A9];
Apim = [Apim1,Apim2,Apim3];

bymid = repmat(by(mid),numpts,1);
by = by - bymid;

A = [Ay;...
    Api;...
    Ai;
    Apim];

b = [by;bpi;bi;pi_m];

% Delete row of zeros. 
A(mid,:) = [];
b(mid,:) = [];

out = A\b;

% Note this is a vecor of 'demeaned EE': y_i - y_mid
y = out(1:numpts);
pi = out(numpts+1:2*numpts);
i = out(2*numpts+1:3*numpts);

Epi = pi'*R.trans_mat(mid,:)';

% Add y_mid back to get y_i
y = y + (pi(mid) - cBET*Epi)/cKAPPA;

out = [y;pi;i];

