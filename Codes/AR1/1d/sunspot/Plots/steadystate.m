function S = steadystate(P,i)

S.pi_targ = P.pi_targ(i);

options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

% Specify Initial Guess
x0 = ones(8,1);
func = @(x) sun_solve(x,P,i);


x_out_sun = fsolve(func,x0,options);


S.c_s    = x_out_sun(1);
S.inf_s  = x_out_sun(2);
S.n_s    = x_out_sun(3);
S.y_s    = S.n_s;
S.w_s    = S.n_s^P.chin*S.c_s^P.chic;
S.r_s    = S.pi_targ/P.beta*((S.inf_s/S.pi_targ)^(P.phi_pi)*(S.y_s/S.y_s)^(P.phi_y));
S.v_s    = x_out_sun(4);

S.c_d    = x_out_sun(5);
S.inf_d  = x_out_sun(6);
S.n_d    = x_out_sun(7);
S.y_d    = S.n_d;
S.w_d    = S.n_d^P.chin*S.c_d^P.chic;
S.r_d    = 1;
S.v_d    = x_out_sun(8);

S.del    = 1;
S.r_zlb  = 1;
S.y = S.y_s;