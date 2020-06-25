%--------------------------------------------------------------------------
% File Name: stylized_rouwenhorst_RAFR_cPHIpi_h.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/ar1/nstate/RAFR
% stylized_rouwenhorst_RAFR_cPHIpi_h
%--------------------------------------------------------------------------


clear all
close all
clc

load('cSIGMAd_max.mat')
%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 10;
cRstar = 1/400;
cRHO    = 0.8;
cSIGMAd_grid = [ 0.15/100,  0.20/100,  cSIGMAd_max.h];
% cSIGMAd_grid = 0.31/100;

%% Housekeeping
nstate = 21;
% Find Middle State
mid = (nstate+1)/2;

% Set Grid Intervals
pi_m_low = -1.5/400;
pi_m_high = 0.5/400;
pts = 1001;
pi_m = linspace(pi_m_low,pi_m_high, pts)';

% Taylor Rule
tr = max(0,cRstar + cPHIpi*pi_m);
% 'Riskless' Fisher Relation
fr = pi_m + cRstar;
% Risk Adjusted Fisher Relation
rafr = zeros(length(pi_m),length(cSIGMAd_grid));


for j = 1:length(cSIGMAd_grid)
    cSIGMAd = cSIGMAd_grid(j);
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd];
    
    %% Define i.i.d Discretizing State-Space
    [R.del_grid,R.trans_mat,R.s]=rouwenhorst(cRHO,cSIGMAd,nstate);
    for k = 1:length(pi_m);
        %% Target Regime
        % Get initial matrix and solution vector
        [A, b, out] = eqmmat(params,R,nstate,pi_m(k));

        y = out(1:nstate);
        pi = out(nstate+1:2*nstate);
        i = out(2*nstate+1:3*nstate);


        % Refine Matrix to account for ZLB
        if sum(i >= 0) == nstate
            converged = 1;        
            rafr(k,j) = cRstar + pi_m(k);
        else
            converged = 0;
        end

        while converged == 0
            [A_up, b_up, out_up] = eqmrefine(params,R,nstate,A,b,pi_m(k));

            i = out_up(2*nstate+1:3*nstate);

            if sum(i >= 0) == nstate
                converged =1;

                y = out_up(1:nstate);
                pi = out_up(nstate+1:2*nstate);
                i = out_up(2*nstate+1:3*nstate);


                Ey = y'*R.trans_mat(mid,:)';
                Epi = pi'*R.trans_mat(mid,:)';

                if sum(i == 0) == nstate
                    rafr(k,j) = cRstar + pi_m(k);
                else
                    rafr(k,j) = cRstar + pi_m(k) + cSIGMA^(-1)*(Ey - y(mid)) + (Epi - pi(mid));
                end

            end

            A = A_up;
            b = b_up;
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
for j = 1:length(cSIGMAd_grid)
    [~,inx] = sort(abs(tr - rafr(:,j)));
    if j == length(cSIGMAd_grid)
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

xlabel('Inflation','FontSize',25)
ylabel('Nominal Interst Rate','FontSize',25)
set(gca,'Xlim',[400*pi_m(1), 400*pi_m(end)],'Ylim',[-1 2],'FontSize',25) 
L = legend([h(1) h(2) h(3)],'\sigma_{\epsilon} = 0.15/100','\sigma_{\epsilon} = 0.20/100','\sigma_{\epsilon} = \sigma_{\epsilon}^{max}');
set(L,'Location','NorthWest','Fontsize',20)

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
% print(fig(1),'-depsc','stylized_rouwenhorst_RAFR_cPHIpi_h.eps');
print(fig(1),'-depsc',strcat(savedir,'stylized_rouwenhorst_RAFR_cPHIpi_h.eps'));

