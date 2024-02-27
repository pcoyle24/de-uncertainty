function C = chebpoly(G,O)

% C = chebpoly(M)
%   Outputs coefficients, powers, and zeros of Chebyshev polynomials
% Inputs:
%     M       :   Structure of Grids
%     O       :   Structure of options
% Output:
%     C.T     :   Chebyshev Polynomial Coefficients
%     C.P     :   Chebyshev Polynomial Powers
%     C.X     :   Chebyshev Polynomial Zeros
%     C.basis :   Chebyshev Basis
    
% Chebyshev polynomial coefficient and power matrix
M = G.griddim;
C.max = max(M)+1;
C.T = zeros(C.max);
C.P = C.T;
C.X = C.T;
for i = 1:C.max
    % Polynomial coefficients
    C.T(i,1:i) = ChebyshevCoeff(i-1);
    for j = 1:i
        % Polynomial powers
        C.P(i,j) = i - j;
        % Zeros of polynomial
        C.X(i,j) = -cos((2*j-1)*pi/(2*i));
        
    end
end

% Calculate Chebyshev basis
C.basis = getbasis1d(O.delbound,G.del,O.n1,C.T,C.P);
% C.basis = getbasis2d(O.sunbound,O.delbound,G.sun,G.del,O.n2,O.n1,C.T,C.P);