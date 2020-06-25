%--------------------------------------------------------------------------
% File Name: RiskAdjustedFisherRelation_syms.m
% Author: Philip Coyle
% Date Created: 01/07/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/3state_iid_shock
% RiskAdjustedFisherRelation_syms
%--------------------------------------------------------------------------

clear all
close all
clc

 %% Model Calibration
syms cBET cSIGMA cKAPPA cPHIpi cRstar cSHOCK p_m pi_m
assume(cBET > 0 & cBET < 1)
assume(cSIGMA > 0)
assume(cKAPPA > 0 & cKAPPA < 1)
assume(cPHIpi > 1)
assume(cRstar > 0)
assume(cSHOCK > 0)
assume(p_m > 0 & p_m < 1)

%% Housekeeping
% Find max inflation that leads to ZLB binding
bound = -cRstar/cPHIpi;
%  Taylor Rule
tr = cRstar + cPHIpi*pi_m;
tr_zlb = 0;
% Standard Fisher Relation
fr = pi_m + cRstar;

tr_dss_pi = solve(tr-fr,pi_m);
dr_dss_pi = solve(tr_zlb-fr,pi_m);

%% Policy Functions from below (Assume middle state interst rate binds)
% Low  High state inflation & expected inflation (as functions of middle state inflation)
pi_l_b = pi_m - cKAPPA*cSIGMA*cSHOCK;
pi_h_b = (pi_m - cKAPPA*cSIGMA*(cRstar - cSHOCK))/(1+cKAPPA*cSIGMA*cPHIpi);
exp_pi_b = (1-p_m)/2*pi_l_b + p_m*pi_m + (1-p_m)/2*pi_h_b;

% Low Middle and High state output & expected output (as functions of middle state inflation)
y_m_b = (pi_m - cBET*exp_pi_b)/cKAPPA;
y_l_b = y_m_b - cSIGMA*cSHOCK;
y_h_b = y_m_b - cSIGMA*(cRstar + cPHIpi*pi_h_b - cSHOCK);
exp_y_b = (1-p_m)/2*y_l_b + p_m*y_m_b + (1-p_m)/2*y_h_b;

% Find pi_m point where high state interst rate binds
lb_cutoff = solve(pi_h_b - bound,pi_m);
dpilb_dcSHOCK = diff(lb_cutoff,cSHOCK);
cSHOCK_lb = solve(lb_cutoff - dr_dss_pi,cSHOCK);

%% Policy Functions from above (Assume middle state interst rate does not bind)
% Low  High state inflation & expected inflation (as functions of middle state inflation)
pi_l_a = pi_m*(1+cKAPPA*cSIGMA*cPHIpi) + cKAPPA*cSIGMA*(cRstar - cSHOCK);
pi_h_a = pi_m + (cKAPPA*cSIGMA*cSHOCK)/(1+cKAPPA*cSIGMA*cPHIpi);
exp_pi_a = (1-p_m)/2*pi_l_a + p_m*pi_m + (1-p_m)/2*pi_h_a;

% Low Middle and High state output & expected output (as functions of middle state inflation)
y_m_a = (pi_m - cBET*exp_pi_a)/cKAPPA;
y_l_a = y_m_a + cSIGMA*(cPHIpi*pi_m + cRstar - cSHOCK);
y_h_a = y_m_a + cSIGMA*(cPHIpi*(pi_m - pi_h_a) + cSHOCK);
exp_y_a = (1-p_m)/2*y_l_a + p_m*y_m_a + (1-p_m)/2*y_h_a;

% Find pi_m point where low state interst rate binds
ub_cutoff = solve(pi_l_a - bound,pi_m);
dpiub_dcSHOCK = diff(ub_cutoff,cSHOCK);
cSHOCK_ub = solve(ub_cutoff - tr_dss_pi,cSHOCK);

%% Define Risk Adjusted Fisher Relationship
rafr_lm = cRstar + pi_m + (exp_pi_b - pi_m) + 1/cSIGMA*(exp_y_b - y_m_b);
rafr_l = cRstar + pi_m + (exp_pi_a - pi_m) + 1/cSIGMA*(exp_y_a - y_m_a);

[rafr_lm_exp,gamma1] = subexpr(collect(simplify(rafr_lm,'steps',1000),pi_m),'gamma1');
[rafr_l_exp,gamma2] = subexpr(collect(simplify(rafr_l,'steps',1000),pi_m),'gamma2');
pretty(rafr_lm_exp);
pretty(rafr_l_exp);

rafr_l_pim = simplify(solve(rafr_l-tr, pi_m),'steps',500);
% rafr_l_pim_zlb = simplify(solve(rafr_l-tr_zlb, pi_m),'steps',500);
% rafr_lm_pim = simplify(solve(rafr_lm-tr, pi_m),'steps',500);
rafr_lm_pim_zlb = simplify(solve(rafr_lm-tr_zlb, pi_m),'steps',500);

dl_dcSHOCK = simplify(diff(rafr_l_pim,cSHOCK),'steps',500);
% dlzlb_dcSHOCK = simplify(diff(rafr_l_pim_zlb,cSHOCK),'steps',500);
% dlm_dcSHOCK = simplify(diff(rafr_lm_pim,cSHOCK),'steps',500);
dlmzlb_dcSHOCK = simplify(diff(rafr_lm_pim_zlb,cSHOCK),'steps',500);


% dl_dcPHIpi = simplify(diff(rafr_l_pim,cPHIpi),'steps',500);
% dlzlb_dcPHIpi = simplify(diff(rafr_l_pim_zlb,cSHOCK),'steps',500);
% dlm_dcPHIpi = simplify(diff(rafr_lm_pim,cSHOCK),'steps',500);
% dlmzlb_dcPHIpi = simplify(diff(rafr_lm_pim_zlb,cPHIpi),'steps',500);

%% Bounds
% Find shock size such that RSS is at bound
rafr_bound = simplify(solve(rafr_l_pim - rafr_lm_pim_zlb, cSHOCK),'steps',500);

% Find cPHIpi bounds
cPHIpi_lb = simplify(solve(cSHOCK_lb - rafr_bound, cPHIpi),'steps',500);
cPHIpi_ub = simplify(solve(cSHOCK_ub - rafr_bound, cPHIpi),'steps',500);


