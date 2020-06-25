function residuals = steady(x_guess,P,i)

C      = x_guess(1);
INF    = x_guess(2);
N      = x_guess(3);
V      = x_guess(4);

PItilde = INF/(P.pi_targ(i)^P.iota*INF^(1-P.iota))^P.alpha;
W = C^P.chic*N^P.chin;
Y = N;
R  = P.pi_targ(i)/P.beta*(INF/P.pi_targ(i))^(P.phi_pi)*(Y/Y).^(P.phi_y);

DELbar = 1;

residuals(1) = C^(-P.chic) - P.beta*DELbar*R*(C^(-P.chic)*INF^(-1));
residuals(2) = N/C^(P.chic)*(P.varphi*(PItilde - 1)*PItilde - (1 - P.theta) - P.theta*(1-P.tau)*W) - P.beta*DELbar*P.varphi*((N/C^(P.chic))*(PItilde - 1)*PItilde);
residuals(3) = N - C - P.varphi/2*(PItilde - 1)^2*N;
if P.chic == 1
    residuals(4) = V - log(C) + N^(1+P.chin)/(1+P.chin) - P.beta*DELbar*V;
else
    residuals(4) = V - C^(1-P.chic)/(1-P.chic) + N^(1+P.chin)/(1+P.chin) - P.beta*DELbar*V;
end


