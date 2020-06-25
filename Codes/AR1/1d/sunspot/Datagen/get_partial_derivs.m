% get_partial_derivs.m
% 
% Analytical Partial Derivatives w.r.t coefficients (for Jacobian Matrix)
% '_dax' corresponds to derivatives w.r.t coefficients related to pf.c.
% '_dbx' corresponds to derivatives w.r.t coefficients related to pf.inf.
% '_dcx' corresponds to derivatives w.r.t coefficients related to pf.v.

dc_dax = @(x) x';
dc_dbx = @(x) 0*x';
dc_dcx = @(x) 0*x';

dpi_dax = @(x) 0*x';
dpi_dbx = @(x) x';
dpi_dcx = @(x) 0*x';

dv_dax = @(x) 0*x';
dv_dbx = @(x) 0*x';
dv_dcx = @(x) x';

dpitilde_dbx = @(x) 1/(S.pi_targ^P.iota*pi_yesterday^(P.iota-1))^P.alpha*dpi_dbx(x);

dn_dax = @(x,pitilde) (1-P.varphi/2*(pitilde-1)^2)^(-1)*dc_dax(x);
dn_dbx = @(x,c,pitilde) -c*(1-(P.varphi/2)*(pitilde-1)^2)^(-2)*-P.varphi*(pitilde-1)*dpitilde_dbx(x);

dy_dax = @(x,pitilde) dn_dax(x,pitilde);
dy_dbx = @(x,c,pitilde) dn_dbx(x,c,pitilde);

dr_dax = @(x,pi,y,pitilde) S.pi_targ/P.beta*(pi/S.pi_targ)^(P.phi_pi)*P.phi_y*(y/S.y)^(P.phi_y-1)*(1/S.y)*dy_dax(x,pitilde);
dr_dbx = @(x,pi,y,c,pitilde) S.pi_targ/P.beta*((P.phi_pi*(pi/S.pi_targ)^(P.phi_pi-1)*(1/S.pi_targ)*dpi_dbx(x)*(y/S.y)^(P.phi_y)) ...
                + ((pi/S.pi_targ)^(P.phi_pi)*P.phi_y*(y/S.y)^(P.phi_y-1)*(1/S.y)*dy_dbx(x,c,pitilde)));
            
drzlb_dax = @(x,pi,y,pitilde) 0*(S.pi_targ/P.beta*(pi/S.pi_targ)^(P.phi_pi)*P.phi_y*(y/S.y)^(P.phi_y-1)*(1/S.y)*dy_dax(x,pitilde));
drzlb_dbx = @(x,pi,y,c,pitilde) 0*(S.pi_targ/P.beta*((P.phi_pi*(pi/S.pi_targ)^(P.phi_pi-1)*(1/S.pi_targ)*dpi_dbx(x)*(y/S.y)^(P.phi_y)) ...
                + ((pi/S.pi_targ)^(P.phi_pi)*P.phi_y*(y/S.y)^(P.phi_y-1)*(1/S.y)*dy_dbx(x,c,pitilde))));            

dw_dax = @(x,c,pitilde,n) P.chin*n^(P.chin-1)*dn_dax(x,pitilde)*c^P.chic + n^P.chin*P.chic*c^(P.chic-1)*dc_dax(x);
dw_dbx = @(x,c,pitilde,n) P.chin*n^(P.chin-1)*dn_dbx(x,c,pitilde)*c^P.chic;