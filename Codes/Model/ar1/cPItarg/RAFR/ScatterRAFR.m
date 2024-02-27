%--------------------------------------------------------------------------
% File Name: stylized_rouwenhorst_RAFR_cPHIpi_m.m
% Author: Philip Coyle
% Date Created: 01/24/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/ar1/nstate/RAFR
% stylized_rouwenhorst_RAFR_cPHIpi_m
%--------------------------------------------------------------------------


clear all
close all
clc


%% Paramaters
cSIGMA = 1.05;
cKAPPA = 0.04;
cPHIpi = 1.75;
cRstar = 0.5/400;
cBET = 1/(1 + cRstar);
cRHO = 0.81; 
cPItarg= 2/400;
cSIGMAd = 0.125/100;


%% Housekeeping
nstate = 21;
% Find Middle State
mid = (nstate+1)/2;

% Set Grid Intervals
pi_m_low = -2.5/400;
pi_m_high = cPItarg + 0.5/400;
pts = 1001;
pi_m = linspace(pi_m_low,pi_m_high, pts)';

% Taylor Rule
tr = max(0,cRstar + cPItarg + cPHIpi*(pi_m - cPItarg));
% 'Riskless' Fisher Relation
fr = pi_m + cRstar;
% Risk Adjusted Fisher Relation
rafr = zeros(length(pi_m),1);


params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cPItarg;cSIGMAd];

%% Define i.i.d Discretizing State-Space
[R.del_grid,R.trans_mat,R.s]=rouwenhorst(cRHO,cSIGMAd,nstate);

% Solve for RAFR
for k = 1:length(pi_m)
    %% Target Regime
    % Get initial matrix and solution vector
    [A, b, out] = eqmmat(params,R,nstate,pi_m(k));

    y = out(1:nstate);
    pi = out(nstate+1:2*nstate);
    i = out(2*nstate+1:3*nstate);


    % Refine Matrix to account for ZLB
    if sum(i >= 0) == nstate
        converged = 1;        
        rafr(k) = cRstar + pi_m(k) - cPItarg;
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


%% Solve for Policy Function
params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd;cPItarg];

% Get initial matrix and solution vector
[A_dr, b_dr, out_dr] = eqmmat_dr(params,R,nstate);

y_dr = out_dr(1:nstate);
pi_dr = out_dr(nstate+1:2*nstate);
i_dr = out_dr(2*nstate+1:3*nstate);

i_real_dr = cRstar + cPItarg + cPHIpi*(pi_dr - cPItarg);

% Refine Matrix to account for ZLB
if sum(i_dr > 0) == sum(i_real_dr > 0)
    converged = 1;
else
    converged = 0;
end
while converged == 0
    [A_up, b_up, out_up] = eqmrefine_dr(params,R,nstate,A_dr,b_dr);

    pi_dr = out_up(nstate+1:2*nstate);
    i_dr = out_up(2*nstate+1:3*nstate);
    i_real_dr = cRstar + cPItarg + cPHIpi*(pi_dr - cPItarg);


    if sum(i_dr > 0) == sum(i_real_dr > 0)
        converged = 1;

        y_dr = 100*out_up(1:nstate);
        pi_dr = 400*out_up(nstate+1:2*nstate);
        i_dr = 400*out_up(2*nstate+1:3*nstate);

        % Check equilibirum existence
        if sum(i_dr<0) > 0
%                 disp('The equilibrium does not exist. Reduce shock size.')
            converged = -1;
        end
    end

    A_dr = A_up;
    b_dr = b_up;
end

% pi_sim = sim_data(params,pi_dr,R,80,1000);
% i_sim = sim_data(params,i_dr,R,80,1000);

        
%% Plotting
fig(1) = figure(1);
grid on
box on
hold on
h(1) = plot(400*pi_m, 400*tr,'Color','k','LineStyle','-','LineWidth',2);
h(2) = plot(400*pi_m, 400*fr,'Color','k','LineStyle','--','LineWidth',2);
% scatter(pi_sim,i_sim,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
for i = 1:nstate
    p = plot(pi_dr(i), i_dr(i));
    set(p,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',100*R.s(i))
end



[~,inx] = sort(abs(tr - rafr(:)));
RSS = inx(1:2);
it = 2;
while abs(RSS(1) - RSS(2)) <= 10
    RSS(2) = inx(it + 1);
    it = it + 1;
end
h(3) = plot(400*pi_m, 400*rafr(:),'Color','k','LineStyle','-.','LineWidth',2);
% p1 = plot(400*pi_m(RSS(1)),400*rafr(RSS(1)));
% p2 = plot(400*pi_m(RSS(2)),400*rafr(RSS(2))); 
% % set([p1 p2],'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8)

xlabel('Inflation','FontSize',25)
ylabel('Nominal Interst Rate','FontSize',25)
set(gca,'XLim',[400*pi_m(1), 400*pi_m(end)],'YLim',[-1, 4],'FontSize',25)
L = legend([h(1) h(2) h(3)],'Taylor Rule','Standard Fisher Relation','Risk Adjusted Fisher Relation');
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
print(fig(1),'-dpdf',strcat(savedir,'scatter_RAFR.pdf'));

