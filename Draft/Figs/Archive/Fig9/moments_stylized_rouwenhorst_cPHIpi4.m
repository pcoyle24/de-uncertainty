%--------------------------------------------------------------------------
% File Name: moments_stylized_rouwenhorst_cPHIpi4.m
% Author: Philip Coyle
% Date Created: 01/18/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/Fig9
% moments_stylized_rouwenhorst_cPHIpi4
%--------------------------------------------------------------------------

clear all
close all
clc


%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 4;
cRstar = 1/400;
cRHO = 0.8; 
cSIGMAd_grid = linspace(0,0.3645/100,51);

numpts = 31;

%% Housekeeping
e_y_tr = zeros(length(cSIGMAd_grid),1);
e_y_dr = zeros(length(cSIGMAd_grid),1);

rss_y_tr = zeros(length(cSIGMAd_grid),1);
rss_y_dr = zeros(length(cSIGMAd_grid),1);

e_pi_tr = zeros(length(cSIGMAd_grid),1);
e_pi_dr = zeros(length(cSIGMAd_grid),1);

rss_pi_tr = zeros(length(cSIGMAd_grid),1);
rss_pi_dr = zeros(length(cSIGMAd_grid),1);

e_i_tr = zeros(length(cSIGMAd_grid),1);
e_i_dr = zeros(length(cSIGMAd_grid),1);

rss_i_tr = zeros(length(cSIGMAd_grid),1);
rss_i_dr = zeros(length(cSIGMAd_grid),1);


e_elb_tr = zeros(length(cSIGMAd_grid),1);
e_elb_dr = zeros(length(cSIGMAd_grid),1);

sim = 1000000;
burn = 1000;

tic;
for i = 1:length(cSIGMAd_grid)
    cSIGMAd = cSIGMAd_grid(i);
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd];


    %% Define Rouwenhorst Discretizing State-Space
    [R.del_grid,R.trans_mat,R.s]=rouwenhorst(cRHO,cSIGMAd,numpts);
    if i == 1;
        rss_inx = (numpts+1)/2;
    else
        rss_inx = find(abs(R.del_grid - 0) < 1e-10);
    end

    %% Target Regime
    % Get initial matrix and solution vector
    [A_tr, b_tr, out_tr] = eqmmat_tr(params,R,numpts);

    y_tr = out_tr(1:numpts);
    pi_tr = out_tr(numpts+1:2*numpts);
    i_tr = out_tr(2*numpts+1:3*numpts);

    % Refine Matrix to account for ZLB
    if sum(i_tr >= 0) == numpts
        converged = 1;
    else
        converged = 0;
    end

    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_tr(params,R,numpts,A_tr,b_tr);

        i_tr = out_up(2*numpts+1:3*numpts);

        if sum(i_tr >= 0) == numpts
            converged =1;

            y_tr = out_up(1:numpts);
            pi_tr = out_up(numpts+1:2*numpts);
            i_tr = out_up(2*numpts+1:3*numpts);

            % Check equilibirum existence
            i_real_tr = cRstar + cPHIpi*pi_tr;
            inx = find(i_real_tr>0);
            if sum(abs(i_tr(inx) - i_real_tr(inx)) < eps) ~= length(inx)
                error('The equilibrium does not exist. Reduce shock size.')
            end
        end

        A_tr = A_up;
        b_tr = b_up;
    end
    
    % Get RSS
    rss_y_tr(i) = y_tr(rss_inx);
    rss_pi_tr(i) = pi_tr(rss_inx);
    rss_i_tr(i) = i_tr(rss_inx);
    
    % Get Expected Values
    e_y_tr(i) = exp_val(params,y_tr,R,sim,burn);
    e_pi_tr(i) = exp_val(params,pi_tr,R,sim,burn);
    [e_i_tr(i), e_elb_tr(i)] = exp_val(params,i_tr,R,sim,burn,1);


    
    %% Deflationary Regime
    % Get initial matrix and solution vector
    [A_dr, b_dr, out_dr] = eqmmat_dr(params,R,numpts);

    y_dr = out_dr(1:numpts);
    pi_dr = out_dr(numpts+1:2*numpts);
    i_dr = out_dr(2*numpts+1:3*numpts);

    i_real_dr = cRstar + cPHIpi*pi_dr;

    % Refine Matrix to account for ZLB
    if sum(i_dr > 0) == sum(i_real_dr > 0)
        converged = 1;
    else
        converged = 0;
    end
    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_dr(params,R,numpts,A_dr,b_dr);

        pi_dr = out_up(numpts+1:2*numpts);
        i_dr = out_up(2*numpts+1:3*numpts);
        i_real_dr = cRstar + cPHIpi*pi_dr;


        if sum(i_dr > 0) == sum(i_real_dr > 0)
            converged =1;

            y_dr = out_up(1:numpts);
            pi_dr = out_up(numpts+1:2*numpts);
            i_dr = out_up(2*numpts+1:3*numpts);

            % Check equilibirum existence
            if sum(i_dr<0) > 0
                error('The equilibrium does not exist. Reduce shock size.');
            end
        end

        A_dr = A_up;
        b_dr = b_up;
    end
    
    % Get RSS
    rss_y_dr(i) = y_dr(rss_inx);
    rss_pi_dr(i) = pi_dr(rss_inx);
    rss_i_dr(i) = i_dr(rss_inx);
    
    % Get Expected Values
    e_y_dr(i) = exp_val(params,y_dr,R,sim,burn);
    e_pi_dr(i) = exp_val(params,pi_dr,R,sim,burn);
    [e_i_dr(i), e_elb_dr(i)] = exp_val(params,i_dr,R,sim,burn,1);
    
    disp(strcat('cSIGMAd = ',num2str(i),'of',num2str(length(cSIGMAd_grid))))
end
toc;

%% Plotting

X_tr = [100*rss_y_tr,400*rss_pi_tr,400*rss_i_tr,100*e_y_tr,400*e_pi_tr,400*e_i_tr,100*e_elb_tr];
X_dr = [100*rss_y_dr,400*rss_pi_dr,400*rss_i_dr,100*e_y_dr,400*e_pi_dr,400*e_i_dr,100*e_elb_dr];
ylims = [0 0.4; -1 0; 0 1;-0.4 0.01; -1.01 0.1; 0 1.2;0 100];

header = {'y_{RSS}','\pi_{RSS}','i_{RSS}','E[y]','E[\pi]','E[i]','Prob[ELB]'};

fig(1) = figure(1);
for i = 1:length(header);
    subplot(3,3,i)
    box on
    grid on
    hold on
    if i == length(header);
        h(1) = plot(cSIGMAd_grid,X_tr(:,i),'Color','k','LineStyle','-','LineWidth',2);
        h(2) = plot(cSIGMAd_grid,X_dr(:,i),'Color','k','LineStyle','--','LineWidth',2);
    else
        plot(cSIGMAd_grid,X_tr(:,i),'Color','k','LineStyle','-','LineWidth',2)
        plot(cSIGMAd_grid,X_dr(:,i),'Color','k','LineStyle','--','LineWidth',2)
    end
    
    xlabel('\sigma_{\epsilon}','FontSize',15)
    title(header{i},'FontSize',15,'FontWeight','Normal')
    set(gca,'Xlim',[cSIGMAd_grid(2) cSIGMAd_grid(end)],'Ylim',ylims(i,:),'FontSize',15)
end
L = legend([h(1) h(2)], 'Target Equilibrium', 'Deflationary Equilibrium');
set(L,'Position',[0.475 0.2 0.1 0.1],'FontSize',15)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','rouwenhorst_moments_cPHIpi4.eps');
