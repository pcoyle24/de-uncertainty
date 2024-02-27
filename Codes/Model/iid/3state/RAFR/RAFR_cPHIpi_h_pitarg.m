%--------------------------------------------------------------------------
% File Name: RAFR_cPHIpi_h_pitarg.m
% Author: Philip Coyle
% Date Created: 01/07/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/iid/3state_iid_shock/RAFR
% RAFR_cPHIpi_h_pitarg
%--------------------------------------------------------------------------

clear all
close all
% clc

%% Model Calibration
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 5;
cRstar = 1/400;
cPItarg = 2/400;
p_m = 0.5; %Probability of being in the middle state 

cSHOCK_lb =((cPItarg + cRstar)*(cPHIpi - 1))/(cKAPPA*cPHIpi*cSIGMA);
cSHOCK_ub =((cPItarg + cRstar)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi*cSIGMA);
cSHOCK_bound = -(2*(cPItarg + cRstar)*(cPHIpi - 1)*(cKAPPA*cPHIpi*cSIGMA + 1))/(cKAPPA*cPHIpi^2*cSIGMA*(cKAPPA*cSIGMA + 1)*(p_m - 1));
cSHOCK_grid = [0.2, 0.275 max([cSHOCK_lb,cSHOCK_ub,cSHOCK_bound])];

%% Housekeeping
% Set Grid Intervals
pi_m_low = -1.5/400;
pi_m_high = 2.5/400;
numgrid = 1001;
pi_m = linspace(pi_m_low,pi_m_high, numgrid)';

% Allocate Space
inx_ub = zeros(length(cSHOCK_grid),1);
pi_ub_zlb = zeros(length(cSHOCK_grid),1);
inx_lb = zeros(length(cSHOCK_grid),1);
pi_lb_zlb = zeros(length(cSHOCK_grid),1);
rafr = zeros(numgrid,length(cSHOCK_grid));


% Find max inflation that leads to ZLB binding
bound = (-cRstar - cPItarg*(1-cPHIpi))/cPHIpi;

inx_m = find(min(abs(pi_m - bound)) == abs(pi_m - bound));
% Find pi_m point where middle state interst rate binds
pi_m_zlb = pi_m(inx_m);

% Taylor Rule
tr = max(0,cRstar + cPItarg + cPHIpi*(pi_m - cPItarg));

% 'Riskless' Fisher Relation
fr = pi_m + cRstar;

for j = 1:length(cSHOCK_grid)
    cSHOCK = cSHOCK_grid(j);
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
    inx_lb(j) = find(min(abs(pi_h_b - bound)) == abs(pi_h_b - bound));
    pi_lb_zlb(j) = pi_m(inx_lb(j));

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
    inx_ub(j) = find(min(abs(pi_l_a - bound)) == abs(pi_l_a - bound));
    pi_ub_zlb(j) = pi_m(inx_ub(j));

    %% Define Risk Adjusted Fisher Relationship   
    for i = 1:numgrid
        if pi_m(i) < pi_lb_zlb(j) || pi_m(i) > pi_ub_zlb(j)
            rafr(i,j) = pi_m(i) + cRstar;
        elseif pi_m(i) <= pi_m_zlb
            rafr(i,j) = cRstar + pi_m(i) + (exp_pi_b(i) - pi_m(i)) + 1/cSIGMA*(exp_y_b(i) - y_m_b(i));
        elseif pi_m(i) <= pi_ub_zlb(j)
            rafr(i,j) = cRstar + pi_m(i) + (exp_pi_a(i) - pi_m(i)) + 1/cSIGMA*(exp_y_a(i) - y_m_a(i));
        end
    end
end
        
%% Plotting
type = {'--','-.',':'};

fig(1) = figure(1);
grid on
box on
hold on
plot(400*pi_m, 400*tr,'Color','k','LineWidth',2);
plot(400*pi_m, 400*fr,'Color','b','LineWidth',2);
for j = 1:length(cSHOCK_grid)
    [~,inx] = sort(abs(tr - rafr(:,j)));
    if j == length(cSHOCK_grid)
        RSS = inx(1);    
    else
        RSS = inx(1:2);
        it = 2;
        while abs(RSS(1) - RSS(2)) <= 10
            RSS(2) = inx(it + 1);
            it = it + 1;
        end
    end
    h(j) = plot(400*pi_m, 400*rafr(:,j),'Color','r','LineStyle',type{j},'LineWidth',2);
    p1 = plot(400*pi_m(RSS(1)),400*rafr(RSS(1),j));
    if j ~= 3
        p2 = plot(400*pi_m(RSS(2)),400*rafr(RSS(2),j)); 
        set([p1 p2],'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8)
    else
        set(p1,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8)
    end
end

xlabel('Inflation','FontSize',16)
ylabel('Nominal Interst Rate','FontSize',16)
set(gca,'Xlim',[400*pi_m(1), 400*pi_m(end)],'Ylim',[-1, 4],'FontSize',25)
L = legend([h(1) h(2) h(3)],'c = 0.20','c = 0.275','c = c_{max}');
set(L,'Location','NorthWest','Fontsize',25)

savedir = cd;
if ispc 
    savedir = fullfile(savedir, '..\..\..');
    savedir = strcat(savedir,'\Final\');
else
    savedir = fullfile(savedir, '../../..');
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','RAFR_cPHIpi_h_pitarg.eps');
print(fig(1),'-depsc',strcat(savedir,'RAFR_cPHIpi_h_pitarg.eps'));

