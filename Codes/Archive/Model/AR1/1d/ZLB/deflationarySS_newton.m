function y = deflationarySS_newton(x)

%Parameters
cBET = 1/1.0075;
cCHIc = 1;
cCHIn = 1;
cTHETA = 10;
cVARPHI = 175;
cPHIpi = 1.5;
cTAU = 0; %1/cTHETA;

pi = x(1);
c  = x(2);
n  = x(3);

% Our root finding problem is given to us as eq 7-8 and 10 in problem set 1 as
% f(x). 

y = [cBET*pi^(-1) - 1;                                                           %Equation 10
    (1-cBET)*(cVARPHI*(pi-1)*pi)-(1-cTHETA)-cTHETA*n^cCHIn*c^cCHIc;              %Equation 7
    n*(1-cVARPHI/2*(pi-1)^2)-c];                                                 %Equation 8
