%--------------------------------------------------------------------------
% File Name: RAFR_derivs.m
% Author: Philip Coyle
% Date Created: 01/07/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/3state_iid_shock
% RAFR_derivs
%--------------------------------------------------------------------------

clear all
close all
clc

cBET = 1./1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cRstar = 1./400;
p_m = 0.5; %Probability of being in the middle state 

cPHIpi_ub = -2/(p_m + cKAPPA*cSIGMA + cKAPPA*cSIGMA*p_m - 1);
cPHIpi_lb = 2/(p_m - cKAPPA*cSIGMA + cKAPPA*cSIGMA*p_m + 1);

% %% cPHIpi in (1, cPHIpi_lb)
% cPHIpi = (1:0.001:cPHIpi_lb)';
% 
% dl_dcSHOCK = (2.*cKAPPA.*cSIGMA.*(p_m - 1))./((cPHIpi.*(p_m - cKAPPA.*cSIGMA + cKAPPA.*cSIGMA.*p_m + 1) - 2).*(p_m + 1)) + (cKAPPA.*cSIGMA.*(p_m - 1))./((p_m + 1).*(cKAPPA.*cPHIpi.*cSIGMA + 1));
% % dlmzlb_dcSHOCK = -(cKAPPA.*cPHIpi.*cSIGMA.*(cKAPPA.*cSIGMA + 1).*(p_m - 1))./(cPHIpi.*p_m - cPHIpi + cKAPPA.*cPHIpi.*cSIGMA + cKAPPA.*cPHIpi.*cSIGMA.*p_m + 2);
% dlmzlb_dcSHOCK = zeros(length(cPHIpi),1);
% 
% fig(1) = figure(1);
% box on
% grid on
% hold on
% h(1) = plot(cPHIpi,dlmzlb_dcSHOCK,'linewidth',2,'color','k');
% h(2) = plot(cPHIpi,dl_dcSHOCK,'linewidth',2,'color','b');
% set(gca,'XLim',[cPHIpi(1) cPHIpi(end)],'YLim',[-0.01 0.5])
% 
% xlabel('\phi_{\pi}','FontSize',16)
% L = legend([h(2) h(1)],'\partialRAFR (L)/\partial\sigma_{r^n}','\partialRAFR (LM ZLB)/\partial\sigma_{r^n}');
% set(L,'Location','NorthWest','Fontsize',16)
% 
% set(fig(1),'PaperOrientation','Landscape');
% set(fig(1),'PaperPosition',[0 0 11 8.5]);
% % print(fig(1),'-depsc','deriv_1cPHIpilb.eps');
% 
% %% cPHIpi in (cPHIpi_lb, cPHIpi_ub)
% cPHIpi = (cPHIpi_lb:0.001:cPHIpi_ub)';
% 
% dl_dcSHOCK = (2.*cKAPPA.*cSIGMA.*(p_m - 1))./((cPHIpi.*(p_m - cKAPPA.*cSIGMA + cKAPPA.*cSIGMA.*p_m + 1) - 2).*(p_m + 1)) + (cKAPPA.*cSIGMA.*(p_m - 1))./((p_m + 1).*(cKAPPA.*cPHIpi.*cSIGMA + 1));
% dlmzlb_dcSHOCK = -(cKAPPA.*cPHIpi.*cSIGMA.*(cKAPPA.*cSIGMA + 1).*(p_m - 1))./(cPHIpi.*p_m - cPHIpi + cKAPPA.*cPHIpi.*cSIGMA + cKAPPA.*cPHIpi.*cSIGMA.*p_m + 2);
% 
% fig(2) = figure(2);
% box on
% grid on
% hold on
% h(1) = plot(cPHIpi,dlmzlb_dcSHOCK,'linewidth',2,'color','k');
% h(2) = plot(cPHIpi,dl_dcSHOCK,'linewidth',2,'color','b');
% set(gca,'XLim',[cPHIpi(1) cPHIpi(end)],'YLim',[-0.5 0.5])
% 
% xlabel('\phi_{\pi}','FontSize',16)
% L = legend([h(2) h(1)],'\partialRAFR (L)/\partial\sigma_{r^n}','\partialRAFR (LM ZLB)/\partial\sigma_{r^n}');
% set(L,'Location','SouthEast','Fontsize',16)
% 
% set(fig(2),'PaperOrientation','Landscape');
% set(fig(2),'PaperPosition',[0 0 11 8.5]);
% % print(fig(2),'-depsc','deriv_cPHIpilbcPHIpiub.eps');
% 
% %% cPHIpi in (cPHIpi_ub, 10)
% cPHIpi = (cPHIpi_ub + 1e-10:0.001:10)';
% 
% % dl_dcSHOCK = (2.*cKAPPA.*cSIGMA.*(p_m - 1))./((cPHIpi.*(p_m - cKAPPA.*cSIGMA + cKAPPA.*cSIGMA.*p_m + 1) - 2).*(p_m + 1)) + (cKAPPA.*cSIGMA.*(p_m - 1))./((p_m + 1).*(cKAPPA.*cPHIpi.*cSIGMA + 1));
% dl_dcSHOCK = zeros(length(cPHIpi),1);
% dlmzlb_dcSHOCK = -(cKAPPA.*cPHIpi.*cSIGMA.*(cKAPPA.*cSIGMA + 1).*(p_m - 1))./(cPHIpi.*p_m - cPHIpi + cKAPPA.*cPHIpi.*cSIGMA + cKAPPA.*cPHIpi.*cSIGMA.*p_m + 2);
% 
% fig(3) = figure(3);
% box on
% grid on
% hold on
% h(1) = plot(cPHIpi,dlmzlb_dcSHOCK,'linewidth',2,'color','k');
% h(2) = plot(cPHIpi,dl_dcSHOCK,'linewidth',2,'color','b');
% set(gca,'XLim',[cPHIpi(1) cPHIpi(end)],'YLim',[-0.5 0.5])
% 
% xlabel('\phi_{\pi}','FontSize',16)
% L = legend([h(2) h(1)],'\partialRAFR (L)/\partial\sigma_{r^n}','\partialRAFR (LM ZLB)/\partial\sigma_{r^n}');
% set(L,'Location','SouthEast','Fontsize',16)
% 
% set(fig(3),'PaperOrientation','Landscape');
% set(fig(3),'PaperPosition',[0 0 11 8.5]);
% % print(fig(3),'-depsc','deriv_cPHIpiub10.eps');

%% cPHIpi in (1, cPHIpi_lb)
cPHIpi = (1:0.001:cPHIpi_lb)';

% dl_dcSHOCK = (2.*cKAPPA.*cSIGMA.*(p_m - 1))./((cPHIpi.*(p_m - cKAPPA.*cSIGMA + cKAPPA.*cSIGMA.*p_m + 1) - 2).*(p_m + 1)) + (cKAPPA.*cSIGMA.*(p_m - 1))./((p_m + 1).*(cKAPPA.*cPHIpi.*cSIGMA + 1));
dl_dcSHOCK = zeros(length(cPHIpi),1);
dlmzlb_dcSHOCK = -(cKAPPA.*cPHIpi.*cSIGMA.*(cKAPPA.*cSIGMA + 1).*(p_m - 1))./(cPHIpi.*p_m - cPHIpi + cKAPPA.*cPHIpi.*cSIGMA + cKAPPA.*cPHIpi.*cSIGMA.*p_m + 2);
% dlmzlb_dcSHOCK = zeros(length(cPHIpi),1);

fig(4) = figure(4);
box on
grid on
hold on
h(1) = plot(cPHIpi,dlmzlb_dcSHOCK,'linewidth',2,'color','k');
h(2) = plot(cPHIpi,dl_dcSHOCK,'linewidth',2,'color','b');
set(gca,'XLim',[cPHIpi(1) cPHIpi(end)],'YLim',[0 0.015])

xlabel('\phi_{\pi}','FontSize',16)
L = legend([h(2) h(1)],'\partialRSS (TR)/\partial\sigma_{r^n}','\partialRSS (DR)/\partial\sigma_{r^n}');
set(L,'Location','NorthWest','Fontsize',16)

set(fig(4),'PaperOrientation','Landscape');
set(fig(4),'PaperPosition',[0 0 11 8.5]);
% print(fig(4),'-depsc','deriv_1cPHIpilb.eps');

%% cPHIpi in (cPHIpi_lb, cPHIpi_ub)
cPHIpi = (cPHIpi_lb:0.001:cPHIpi_ub)';

dl_dcSHOCK = (2.*cKAPPA.*cSIGMA.*(p_m - 1))./((cPHIpi.*(p_m - cKAPPA.*cSIGMA + cKAPPA.*cSIGMA.*p_m + 1) - 2).*(p_m + 1)) + (cKAPPA.*cSIGMA.*(p_m - 1))./((p_m + 1).*(cKAPPA.*cPHIpi.*cSIGMA + 1));
dlmzlb_dcSHOCK = -(cKAPPA.*cPHIpi.*cSIGMA.*(cKAPPA.*cSIGMA + 1).*(p_m - 1))./(cPHIpi.*p_m - cPHIpi + cKAPPA.*cPHIpi.*cSIGMA + cKAPPA.*cPHIpi.*cSIGMA.*p_m + 2);

fig(5) = figure(5);
box on
grid on
hold on
h(1) = plot(cPHIpi,dlmzlb_dcSHOCK,'linewidth',2,'color','k');
h(2) = plot(cPHIpi,dl_dcSHOCK,'linewidth',2,'color','b');
set(gca,'XLim',[cPHIpi(1) cPHIpi(end)],'YLim',[-0.5 0.5])

xlabel('\phi_{\pi}','FontSize',16)
L = legend([h(2) h(1)],'\partialRSS (TR)/\partial\sigma_{r^n}','\partialRSS (DR)/\partial\sigma_{r^n}');
set(L,'Location','SouthEast','Fontsize',16)

set(fig(5),'PaperOrientation','Landscape');
set(fig(5),'PaperPosition',[0 0 11 8.5]);
% print(fig(5),'-depsc','deriv_cPHIpilbcPHIpiub.eps');

%% cPHIpi in (cPHIpi_ub, 10)
cPHIpi = (cPHIpi_ub + 1e-10:0.001:10)';

dl_dcSHOCK = (2.*cKAPPA.*cSIGMA.*(p_m - 1))./((cPHIpi.*(p_m - cKAPPA.*cSIGMA + cKAPPA.*cSIGMA.*p_m + 1) - 2).*(p_m + 1)) + (cKAPPA.*cSIGMA.*(p_m - 1))./((p_m + 1).*(cKAPPA.*cPHIpi.*cSIGMA + 1));
% dl_dcSHOCK = zeros(length(cPHIpi),1);
% dlmzlb_dcSHOCK = -(cKAPPA.*cPHIpi.*cSIGMA.*(cKAPPA.*cSIGMA + 1).*(p_m - 1))./(cPHIpi.*p_m - cPHIpi + cKAPPA.*cPHIpi.*cSIGMA + cKAPPA.*cPHIpi.*cSIGMA.*p_m + 2);
dlmzlb_dcSHOCK = zeros(length(cPHIpi),1);

fig(6) = figure(6);
box on
grid on
hold on
h(1) = plot(cPHIpi,dlmzlb_dcSHOCK,'linewidth',2,'color','k');
h(2) = plot(cPHIpi,dl_dcSHOCK,'linewidth',2,'color','b');
set(gca,'XLim',[cPHIpi(1) cPHIpi(end)],'YLim',[-0.015 0])

xlabel('\phi_{\pi}','FontSize',16)
L = legend([h(2) h(1)],'\partialRSS (TR)/\partial\sigma_{r^n}','\partialRSS (DR)/\partial\sigma_{r^n}');
set(L,'Location','SouthEast','Fontsize',16)

set(fig(6),'PaperOrientation','Landscape');
set(fig(6),'PaperPosition',[0 0 11 8.5]);
% print(fig(6),'-depsc','deriv_cPHIpiub10.eps');

