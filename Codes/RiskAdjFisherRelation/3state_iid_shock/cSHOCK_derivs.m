clear all
close all
clc

 %% Model Calibration
syms cBET cSIGMA cKAPPA cPHIpi cRstar cSHOCK p_m pi_m cPItarg
assume(cBET > 0 & cBET < 1)
assume(cSIGMA > 0)
assume(cPHIpi > 1)
assume(cKAPPA > 0 & cKAPPA < 1)
assume(cRstar > 0)
assume(cSHOCK > 0)
assume(p_m > 0 & p_m < 1)
assume(cPItarg > 0);


%% Get Derivatives w.r.t. cPHIpi
cSHOCK_lb =((cPItarg + cRstar)*(cPHIpi - 1))/(cKAPPA*cPHIpi*cSIGMA);
dcSHOCKlb_dcPHIpi = diff(cSHOCK_lb,cPHIpi);

cSHOCK_ub =((cPItarg + cRstar)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi*cSIGMA);
dcSHOCKub_dcPHIpi = diff(cSHOCK_ub,cPHIpi);

cSHOCK_bound = -(2*(cRstar + cPItarg)*(cPHIpi - 1)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi^2*cSIGMA*(cKAPPA*cSIGMA + 1)*(p_m - 1));
dcSHOCKbound_dcPHIpi = simplify(diff(cSHOCK_bound,cPHIpi),'steps',500);


%% Let cPHIpi = cPHIpi_lb
cPHIpi_lb = 2/(p_m - cKAPPA*cSIGMA + cKAPPA*cSIGMA*p_m + 1);

cSHOCK_lb = simplify(subs(cSHOCK_lb,cPHIpi,cPHIpi_lb),'steps',500);
cSHOCK_ub = simplify(subs(cSHOCK_ub,cPHIpi,cPHIpi_lb),'steps',500);
cSHOCK_bound = simplify(subs(cSHOCK_bound,cPHIpi,cPHIpi_lb),'steps',500);

% Notice that cSHOCK_ub = cSHOCK_bound. We look at the derivative at this
% point to determine, for a small negative change in epsilon, which has the
% larget shock size.
dcSHOCKlb_dcPHIpi = simplify(subs(dcSHOCKlb_dcPHIpi,cPHIpi,cPHIpi_lb),'steps',500);
dcSHOCKub_dcPHIpi = simplify(subs(dcSHOCKub_dcPHIpi,cPHIpi,cPHIpi_lb),'steps',500);
dcSHOCKbound_dcPHIpi = simplify(subs(dcSHOCKbound_dcPHIpi,cPHIpi,cPHIpi_lb),'steps',500);


%% Let cPHIpi = cPHIpi_ub

cPHIpi_ub = -2/(p_m + cKAPPA*cSIGMA + cKAPPA*cSIGMA*p_m - 1);

cSHOCK_lb = simplify(subs(cSHOCK_lb,cPHIpi,cPHIpi_lb),'steps',500);
cSHOCK_ub = simplify(subs(cSHOCK_ub,cPHIpi,cPHIpi_lb),'steps',500);
cSHOCK_bound = simplify(subs(cSHOCK_bound,cPHIpi,cPHIpi_lb),'steps',500);


% Notice that cSHOCK_lb = cSHOCK_bound. We look at the derivative at this
% point to determine, for a small negative change in epsilon, which has the
% larget shock size.
dcSHOCKlb_dcPHIpi = simplify(subs(dcSHOCKlb_dcPHIpi,cPHIpi,cPHIpi_ub,'steps',500);
dcSHOCKub_dcPHIpi = simplify(subs(dcSHOCKub_dcPHIpi,cPHIpi,cPHIpi_ub,'steps',500);
dcSHOCKbound_dcPHIpi = simplify(subs(dcSHOCKbound_dcPHIpi,cPHIpi,cPHIpi_ub),'steps',500);