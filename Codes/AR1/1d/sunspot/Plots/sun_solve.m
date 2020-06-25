function residuals_sun = sun_solve(x_guess_sun,P,i)

Cs      = x_guess_sun(1);
INFs    = x_guess_sun(2);
Ns      = x_guess_sun(3);
Vs      = x_guess_sun(4); 
Ys      = Ns;

PItildes = INFs/(P.pi_targ(i)^P.iota*INFs^(1-P.iota))^P.alpha;
Ws = Cs^P.chic*Ns^P.chin;

Cd      = x_guess_sun(5);
INFd    = x_guess_sun(6);
Nd      = x_guess_sun(7);
Vd      = x_guess_sun(8); 
Yd      = Nd;

PItilded = INFd/(P.pi_targ(i)^P.iota*INFd^(1-P.iota))^P.alpha;
Wd = Cd^P.chic*Nd^P.chin;

Rs  = P.pi_targ(i)/P.beta*(INFs/P.pi_targ(i))^(P.phi_pi)*(Ys/Ys).^(P.phi_y);
Rd  = 1;

DELbar = 1;

residuals_sun(1) = Cs^(-P.chic) - P.beta*DELbar*Rs*(P.Ps*(Cs^(-P.chic)*INFs^(-1)) + (1-P.Ps)*(Cd^(-P.chic)*INFd^(-1)));
residuals_sun(2) = Ns/Cs^(P.chic)*(P.varphi*(PItildes - 1)*PItildes - (1 - P.theta) - P.theta*(1-P.tau)*Ws) - (P.beta*DELbar*P.varphi*(P.Ps*((Ns/Cs^(P.chic))*(PItildes - 1)*PItildes) + (1-P.Ps)*((Nd/Cd^(P.chic))*(PItilded - 1)*PItilded)));
residuals_sun(3) = Ns - Cs - P.varphi/2*(PItildes - 1)^2*Ns;
if P.chic == 1
    residuals_sun(4) = Vs - (log(Cs) - Ns^(1+P.chin)/(1+P.chin) + P.beta*DELbar*(P.Ps*Vs + (1-P.Ps)*Vd));
else
    residuals_sun(4) = Vs - (Cs^(1-P.chic)/(1-P.chic) - Ns^(1+P.chin)/(1+P.chin) + P.beta*DELbar*(P.Ps*Vs + (1-P.Ps)*Vd));
end

residuals_sun(5) = Cd^(-P.chic) - P.beta*DELbar*Rd*(P.Pd*(Cd^(-P.chic)*INFd^(-1)) + (1-P.Pd)*(Cs^(-P.chic)*INFs^(-1)));
residuals_sun(6) = Nd/Cd^(P.chic)*(P.varphi*(PItilded - 1)*PItilded - (1 - P.theta) - P.theta*(1-P.tau)*Wd) - (P.beta*DELbar*P.varphi*(P.Pd*((Nd/Cd^(P.chic))*(PItilded - 1)*PItilded) + (1-P.Pd)*((Ns/Cs^(P.chic))*(PItildes - 1)*PItildes)));
residuals_sun(7) = Nd - Cd - P.varphi/2*(PItilded - 1)^2*Nd;
if P.chic == 1
    residuals_sun(8) = Vd - (log(Cd) - Nd^(1+P.chin)/(1+P.chin) + P.beta*DELbar*(P.Pd*Vd + (1-P.Pd)*Vs));
else
    residuals_sun(8) = Vs - (Cd^(1-P.chic)/(1-P.chic) - Nd^(1+P.chin)/(1+P.chin) + P.beta*DELbar*(P.Pd*Vd + (1-P.Pd)*Vs));
end

