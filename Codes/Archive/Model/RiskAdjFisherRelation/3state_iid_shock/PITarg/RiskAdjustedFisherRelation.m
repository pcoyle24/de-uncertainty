%--------------------------------------------------------------------------
% File Name: RiskAdjustedFisherRelation.m
% Author: Philip Coyle
% Date Created: 01/22/2019
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/3state_iid_shock/PITarg
% RiskAdjustedFisherRelation
%--------------------------------------------------------------------------

clear all
close all
clc

%% Model Calibration
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 2;
cRstar = 1/400;
cSHOCK = 0.15;
cPItarg = 2/400;
p_m = 0.5; %Probability of being in the middle state 

%% Housekeeping
% Set Grid Intervals
pi_m_low = -1.5/400;
pi_m_high = 2.5/400;
numgrid = 1001;
pi_m = linspace(pi_m_low,pi_m_high, numgrid)';

% Find max inflation that leads to ZLB binding
bound = (-cRstar - cPItarg*(1-cPHIpi))/cPHIpi;

inx_m = find(min(abs(pi_m - bound)) == abs(pi_m - bound));
% Find pi_m point where middle state interst rate binds
if length(inx_m ) > 1
    pi_m_zlb = min(pi_m(inx_m));
else
    pi_m_zlb = pi_m(inx_m);
end

% Taylor Rule
tr = max(0,cRstar + cPItarg + cPHIpi*(pi_m - cPItarg));

% 'Riskless' Fisher Relation
fr = pi_m + cRstar;

%% Policy Functions from below (Assume middle state interst rate binds)
% Low  High state inflation & expected inflation (as functions of middle state inflation)
pi_l_b = pi_m - cKAPPA*cSIGMA*cSHOCK;
pi_h_b = (pi_m - cKAPPA*cSIGMA*(cRstar + cPItarg*(1-cPHIpi)- cSHOCK))/(1+cKAPPA*cSIGMA*cPHIpi);
exp_pi_b = (1-p_m)/2*pi_l_b + p_m*pi_m + (1-p_m)/2*pi_h_b;

% Low Middle and High state output & expected output (as functions of middle state inflation)
y_m_b = (pi_m - cBET*exp_pi_b)/cKAPPA;
y_l_b = y_m_b - cSIGMA*cSHOCK;
y_h_b = y_m_b - cSIGMA*(cRstar + cPItarg + cPHIpi*(pi_h_b - cPItarg) - cSHOCK);
exp_y_b = (1-p_m)/2*y_l_b + p_m*y_m_b + (1-p_m)/2*y_h_b;

% Find pi_m point where high state interst rate binds
inx_lb = find(min(abs(pi_h_b - bound)) == abs(pi_h_b - bound));
if length(inx_m ) > 1
    pi_lb_zlb = min(pi_m(inx_lb));
else
    pi_lb_zlb = pi_m(inx_lb);
end

%% Policy Functions from above (Assume middle state interst rate does not bind)
% Low  High state inflation & expected inflation (as functions of middle state inflation)
pi_l_a = pi_m*(1+cKAPPA*cSIGMA*cPHIpi) + cKAPPA*cSIGMA*(-cPHIpi*cPItarg + cRstar + cPItarg - cSHOCK);
pi_h_a = pi_m + (cKAPPA*cSIGMA*cSHOCK)/(1+cKAPPA*cSIGMA*cPHIpi);
exp_pi_a = (1-p_m)/2*pi_l_a + p_m*pi_m + (1-p_m)/2*pi_h_a;

% Low Middle and High state output & expected output (as functions of middle state inflation)
y_m_a = (pi_m - cBET*exp_pi_a)/cKAPPA;
y_l_a = y_m_a + cSIGMA*(cPHIpi*(pi_m - cPItarg) + cRstar + cPItarg - cSHOCK);
y_h_a = y_m_a + cSIGMA*(cPHIpi*(pi_m - pi_h_a) + cSHOCK);
exp_y_a = (1-p_m)/2*y_l_a + p_m*y_m_a + (1-p_m)/2*y_h_a;

% Find pi_m point where low state interst rate binds
inx_ub = find(min(abs(pi_l_a - bound)) == abs(pi_l_a - bound));
if length(inx_ub) > 1
    pi_ub_zlb = min(pi_m(inx_ub));
else
    pi_ub_zlb = pi_m(inx_ub);
end

%% Define Risk Adjusted Fisher Relationship
rafr = zeros(numgrid,1);

for i = 1:numgrid
    if pi_m(i) < pi_lb_zlb || pi_m(i) > pi_ub_zlb
        rafr(i) = pi_m(i) + cRstar;
    elseif pi_m(i) <= pi_m_zlb
        rafr(i) = cRstar + pi_m(i) + (exp_pi_b(i) - pi_m(i)) + 1/cSIGMA*(exp_y_b(i) - y_m_b(i));
    elseif pi_m(i) <= pi_ub_zlb
        rafr(i) = cRstar + pi_m(i) + (exp_pi_a(i) - pi_m(i)) + 1/cSIGMA*(exp_y_a(i) - y_m_a(i));
    end
end

[~,inx] = sort(abs(tr - rafr));
RSS = inx(1:2);
%% Plotting
fig(1) = figure(1);
grid on
box on
hold on
h(1) = plot(400*pi_m, 400*tr,'Color','k','LineWidth',2);
h(2) = plot(400*pi_m, 400*fr,'Color','b','LineWidth',2);
h(3) = plot(400*pi_m, 400*rafr,'Color','r','LineStyle',':','LineWidth',2);
h1 = plot(400*pi_m(RSS(1)),400*rafr(RSS(1)));
h2 = plot(400*pi_m(RSS(2)),400*rafr(RSS(2)));
set([h1 h2],'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',8)

xlabel('Inflation','FontSize',16)
ylabel('Nominal Interst Rate','FontSize',16)
set(gca,'FontSize',16)
L = legend([h(1) h(2) h(3)],'Taylor Rule','Standard Fisher Relation','Risk-Adjusted Fisher Relation');
set(L,'Location','NorthWest','Fontsize',15)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','RAFR.eps');

