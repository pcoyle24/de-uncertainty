function tk = ChebyshevPoly(n)

% ChebyshevPoly computes the Chebychev polynomial coefficients
% written by David Terr, Raytheon, 5-10-04

% Inputs:
%   n      : Order of the Chebychev polynomial
% Outputs:
%  tk      : Coefficient vector of the Chebychev polynomial

if n==0 
    tk = 1;
elseif n==1
    tk = [1 0]';
else    
    tkm2 = zeros(n+1,1);
    tkm2(n+1) = 1;
    tkm1 = zeros(n+1,1);
    tkm1(n) = 1;
    for k=2:n        
        tk = zeros(n+1,1);
        for e=n-k+1:2:n
            tk(e) = 2*tkm1(e+1) - tkm2(e);
        end        
        if mod(k,2)==0
            tk(n+1) = (-1)^(k/2);
        end        
        if k<n
            tkm2 = tkm1;
            tkm1 = tk;
        end        
    end    
end