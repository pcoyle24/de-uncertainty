function x = quadroot(a,b,c)

% Compute the roots of a quadratic function of the form 
% a*x^2+b*x+c = 0


x1 = (-b + (b.^2 - 4.*a.*c).^(1/2))./(2.*a);
x2 = (-b - (b.^2 - 4.*a.*c).^(1/2))./(2.*a);

x = [x1,x2];
