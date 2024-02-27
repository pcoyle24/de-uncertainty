function residuals_param = param_solve_cSIGMAd(x_guess_param,params,nstate,rss_inx)

cSIGMAd = x_guess_param;
cBET    = params(1);
cSIGMA  = params(2);
cKAPPA  = params(3);
cPHIpi  = params(4);
cRstar  = params(5);

params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd];


%% Define i.i.d Discretizing State-Space
[R.del_grid,R.s]=iid(cSIGMAd,nstate);
    
%% Target Regime
% Get initial matrix and solution vector
[A_tr, b_tr, out_tr] = eqmmat_tr(params,R,nstate);

y_tr = out_tr(1:nstate);
pi_tr = out_tr(nstate+1:2*nstate);
i_tr = out_tr(2*nstate+1:3*nstate);

% Refine Matrix to account for ZLB
if sum(i_tr >= 0) == nstate
    converged = 1;
else
    converged = 0;
end

while converged == 0
    [A_up, b_up, out_up] = eqmrefine_tr(params,R,nstate,A_tr,b_tr);

    i_tr = out_up(2*nstate+1:3*nstate);

    if sum(i_tr >= 0) == nstate
        converged =1;

        y_tr = out_up(1:nstate);
        pi_tr = out_up(nstate+1:2*nstate);
        i_tr = out_up(2*nstate+1:3*nstate);

%         % Check equilibirum existence
%         i_real_tr = cRstar + cPHIpi*pi_tr;
%         inx = find(i_real_tr>0);
%         if sum(abs(i_tr(inx) - i_real_tr(inx)) < eps) ~= length(inx)
%             error('The equilibrium does not exist. Reduce shock size.')
%         end
    end

    A_tr = A_up;
    b_tr = b_up;
end

% Get RSS
rss_pi_tr = pi_tr(rss_inx);

%% Deflationary Regime
% Get initial matrix and solution vector
[A_dr, b_dr, out_dr] = eqmmat_dr(params,R,nstate);

y_dr = out_dr(1:nstate);
pi_dr = out_dr(nstate+1:2*nstate);
i_dr = out_dr(2*nstate+1:3*nstate);

i_real_dr = cRstar + cPHIpi*pi_dr;

% Refine Matrix to account for ZLB
if sum(i_dr > 0) == sum(i_real_dr > 0)
    converged = 1;
else
    converged = 0;
end
while converged == 0
    [A_up, b_up, out_up] = eqmrefine_dr(params,R,nstate,A_dr,b_dr);

    pi_dr = out_up(nstate+1:2*nstate);
    i_dr = out_up(2*nstate+1:3*nstate);
    i_real_dr = cRstar + cPHIpi*pi_dr;


    if sum(i_dr > 0) == sum(i_real_dr > 0)
        converged =1;

        y_dr = out_up(1:nstate);
        pi_dr = out_up(nstate+1:2*nstate);
        i_dr = out_up(2*nstate+1:3*nstate);

    end

    A_dr = A_up;
    b_dr = b_up;
end

% Get RSS
rss_pi_dr = pi_dr(rss_inx);

%% Residuals
residuals_param(1) = rss_pi_tr - rss_pi_dr;
