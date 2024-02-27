%--------------------------------------------------------------------------
% File Name: stylized_iid.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/AR1/Rouwenhorst
% stylized_iid
%--------------------------------------------------------------------------

clear all
close all
clc

%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 2;
cRstar = 1/400;
cSIGMAd = 0.05;

params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd];

numpts = 3;
% Find Middle State
mid = (numpts+1)/2;

pi_m = linspace(-1.5/400,0.5/400,101);
% pi_m = 0.5/400;

% Taylor Rule
tr = max(0,cRstar + cPHIpi*pi_m);
% 'Riskless' Fisher Relation
fr = pi_m + cRstar;
% Risk Adjusted Fisher Relation
rafr = zeros(length(pi_m),1);

%% Define i.i.d Discretizing State-Space
[R.del_grid,R.s]=iid(cSIGMAd,1,numpts);

for k = 1:length(pi_m);
    %% Target Regime
    % Get initial matrix and solution vector
    [A, b, out] = eqmmat(params,R,numpts,pi_m(k));

    y = out(1:numpts);
    pi = out(numpts+1:2*numpts);
    i = out(2*numpts+1:3*numpts);


    % Refine Matrix to account for ZLB
    if sum(i >= 0) == numpts
        converged = 1;        
        rafr(k) = cRstar + pi_m(k);
    else
        converged = 0;
    end

    while converged == 0
        [A_up, b_up, out_up] = eqmrefine(params,R,numpts,A,b,pi_m(k));

        i = out_up(2*numpts+1:3*numpts);

        if sum(i >= 0) == numpts
            converged =1;

            y = out_up(1:numpts);
            pi = out_up(numpts+1:2*numpts);
            i = out_up(2*numpts+1:3*numpts);

            
            Ey = y'*R.s;
            Epi = pi'*R.s;
            
            if sum(i == 0) == numpts
                rafr(k) = cRstar + pi_m(k);
            else
                rafr(k) = cRstar + pi_m(k) + cSIGMA^(-1)*(Ey - y(mid)) + (Epi - pi(mid));
%                 if cRstar + pi_m(k) + cSIGMA^(-1)*(Ey - y(mid)) + (Epi - pi(mid)) > cRstar + pi_m(k)
%                     error('See here')
%                 end
            end
            
        end

        A = A_up;
        b = b_up;
    end
end

%% Plotting
fig(1) = figure(1);
grid on
box on
hold on
h(1) = plot(400*pi_m, 400*tr,'Color','k','LineWidth',2);
h(2) = plot(400*pi_m, 400*fr,'Color','b','LineWidth',2);
h(3) = plot(400*pi_m, 400*rafr,'Color','r','LineStyle',':','LineWidth',2);
% h1 = plot(400*pi_m(RSS(1)),400*rafr(RSS(1)));
% h2 = plot(400*pi_m(RSS(2)),400*rafr(RSS(2)));
% set([h1 h2],'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',8)

xlabel('Inflation','FontSize',16)
ylabel('Nominal Interst Rate','FontSize',16)
set(gca,'FontSize',16)
L = legend([h(1) h(2) h(3)],'Taylor Rule','Standard Fisher Relation','Risk-Adjusted Fisher Relation');
set(L,'Location','NorthWest','Fontsize',15)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','RAFR.eps');


