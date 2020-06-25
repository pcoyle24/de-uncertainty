%--------------------------------------------------------------------------
% File Name: stylized_rouwenhorst_RAFR.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/ar1/nstate/RAFR
% stylized_rouwenhorst_RAFR
%--------------------------------------------------------------------------

clear all
close all
clc

%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 3;
cRstar = 1/400;
cSIGMAd = 0.19/100;
cRHO    = 0.8;

params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd];

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
rafr = zeros(length(pi_m),1);

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
        rafr(k) = cRstar + pi_m(k);
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
                rafr(k) = cRstar + pi_m(k);
            else
                rafr(k) = cRstar + pi_m(k) + cSIGMA^(-1)*(Ey - y(mid)) + (Epi - pi(mid));
            end
            
        end

        A = A_up;
        b = b_up;
    end
end

%% Get DSS and RSS
[~,inx] = sort(abs(tr - fr));
DSS = inx(1:2);
it = 2;
while abs(DSS(1) - DSS(2)) <= 10
    DSS(2) = inx(it + 1);
    it = it + 1;
end



[~,inx] = sort(abs(tr - rafr));
RSS = inx(1:2);
it = 2;
while abs(RSS(1) - RSS(2)) <= 10
    RSS(2) = inx(it + 1);
    it = it + 1;
end

%% Plotting
fig(1) = figure(1);
grid off
box off
hold on

% Plot Lines
h(1) = plot(400*pi_m, 400*tr,'Color','k','LineWidth',2);
h(2) = plot(400*pi_m, 400*fr,'Color','b','LineWidth',2);
h(3) = plot(400*pi_m, 400*rafr,'Color','r','LineStyle','-','LineWidth',2);

% Plot DSS
pdss1 = plot(400*pi_m(DSS(1)),400*tr(DSS(1)));
pdss2 = plot(400*pi_m(DSS(2)),400*tr(DSS(2))); 
set(pdss1,'Marker','d','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12)
set(pdss2,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12)

% Plot RSS
prss1 = plot(400*pi_m(RSS(1)),400*rafr(RSS(1)));
prss2 = plot(400*pi_m(RSS(2)),400*rafr(RSS(2))); 
set(prss1,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12)
set(prss2,'Marker','d','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12)

% Clean up x and y axis
set(gca,'XLim',[400*pi_m(1), 400*pi_m(end)],'XTick',[-400*cRstar 0],'XTickLabel',{'-r*','\pi*' },'YLim',[-1 3],'YTick',[0],'YTickLabel',{'0'},'FontSize',25)
inf = text(400*pi_m(end),-1.2,'\bf \pi','Interpreter','tex');
R = text(-1.6,2.9,'\bf R','Interpreter','tex');
set(inf,'FontWeight','bold','HorizontalAlignment','center','FontSize',30)
set(R,'FontWeight','bold','HorizontalAlignment','center','FontSize',30)
annotation('arrow',[0.13 .92],[0.11 0.11])
annotation('arrow',[0.13 0.13],[0.11 0.94])

% Make legend
L = legend([h(1) h(2) h(3) pdss1 pdss2 prss1 prss2],'Taylor Rule','Standard Fisher Relation','Risk-Adjusted Fisher Relation','DSS (TE)','DSS (DE)','RSS (TE)','RSS (DE)');
set(L,'Location','NorthWest','Fontsize',20)
set([pdss1 pdss2 prss1 prss2],'linestyle','none')
    legend('boxoff')


%% Saving Figure
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
% print(fig(1),'-depsc','stylized_rouwenhorst_RAFR.eps');
print(fig(1),'-depsc',strcat(savedir,'RAFR.eps'));



