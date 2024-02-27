clear all
close all
clc

 %% Model Calibration
syms cBET cSIGMA cKAPPA cPHIpi cRstar cSHOCK p_m pi_m 
assume(cBET > 0 & cBET < 1)
assume(cSIGMA > 0)
assume(cPHIpi > 1)
assume(cKAPPA > 0 & cKAPPA < 1)
assume(cRstar > 0)
assume(cSHOCK > 0)
assume(p_m > 0 & p_m < 1)


%% Get Derivatives w.r.t. cPHIpi
cSHOCK_lb =(cRstar - cRstar/cPHIpi)/(cKAPPA*cSIGMA);
dcSHOCKlb_dcPHIpi = diff(cSHOCK_lb,cPHIpi);

cSHOCK_ub =(cRstar/cPHIpi + cKAPPA*cRstar*cSIGMA)/(cKAPPA*cSIGMA);
dcSHOCKub_dcPHIpi = diff(cSHOCK_ub,cPHIpi);

cSHOCK_bound = -(2*cRstar*(cPHIpi - 1)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi^2*cSIGMA*(cKAPPA*cSIGMA + 1)*(p_m - 1));
dcSHOCKbound_dcPHIpi = simplify(diff(cSHOCK_bound,cPHIpi),'steps',500);


%% Let cPHIpi = cPHIpi_lb
cPHIpi_lb = 2/(p_m - cKAPPA*cSIGMA + cKAPPA*cSIGMA*p_m + 1);

cSHOCK_lb =(cRstar - cRstar/cPHIpi_lb)/(cKAPPA*cSIGMA);
cSHOCK_ub =(cRstar/cPHIpi_lb + cKAPPA*cRstar*cSIGMA)/(cKAPPA*cSIGMA);
cSHOCK_bound = -(2*cRstar*(cPHIpi_lb - 1)*(cKAPPA*cPHIpi_lb*cSIGMA + 1))/(cKAPPA*cPHIpi_lb^2*cSIGMA*(cKAPPA*cSIGMA + 1)*(p_m - 1));

% Notice that cSHOCK_ub > cSHOCK_lb
cSHOCK_lb = simplify(cSHOCK_lb,'steps',1000);
cSHOCK_ub = simplify(cSHOCK_ub,'steps',500);
cSHOCK_bound = simplify(cSHOCK_bound,'steps',500);

% Notice that cSHOCK_ub = cSHOCK_bound. We look at the derivative at this
% point to determine, for a small negative change in epsilon, which has the
% larget shock size.
dcSHOCKlb_dcPHIpi = subs(dcSHOCKlb_dcPHIpi,cPHIpi,cPHIpi_lb);
dcSHOCKub_dcPHIpi = subs(dcSHOCKub_dcPHIpi,cPHIpi,cPHIpi_lb);
dcSHOCKbound_dcPHIpi = simplify(subs(dcSHOCKbound_dcPHIpi,cPHIpi,cPHIpi_lb),'steps',500);


%% Let cPHIpi = cPHIpi_ub

cPHIpi_ub = -2/(p_m + cKAPPA*cSIGMA + cKAPPA*cSIGMA*p_m - 1);

cSHOCK_lb =(cRstar - cRstar/cPHIpi_ub)/(cKAPPA*cSIGMA);
cSHOCK_ub =(cRstar/cPHIpi_ub + cKAPPA*cRstar*cSIGMA)/(cKAPPA*cSIGMA);
cSHOCK_bound = -(2*cRstar*(cPHIpi_ub - 1)*(cKAPPA*cPHIpi_ub*cSIGMA + 1))/(cKAPPA*cPHIpi_ub^2*cSIGMA*(cKAPPA*cSIGMA + 1)*(p_m - 1));

% Notice that cSHOCK_lb > cSHOCK_ub
cSHOCK_lb = simplify(cSHOCK_lb,'steps',1000);
cSHOCK_ub = simplify(cSHOCK_ub,'steps',500);
cSHOCK_bound = simplify(cSHOCK_bound,'steps',500);

% Notice that cSHOCK_lb = cSHOCK_bound. We look at the derivative at this
% point to determine, for a small negative change in epsilon, which has the
% larget shock size.
dcSHOCKlb_dcPHIpi = subs(dcSHOCKlb_dcPHIpi,cPHIpi,cPHIpi_ub);
dcSHOCKub_dcPHIpi = subs(dcSHOCKub_dcPHIpi,cPHIpi,cPHIpi_ub);
dcSHOCKbound_dcPHIpi = simplify(subs(dcSHOCKbound_dcPHIpi,cPHIpi,cPHIpi_ub),'steps',500);