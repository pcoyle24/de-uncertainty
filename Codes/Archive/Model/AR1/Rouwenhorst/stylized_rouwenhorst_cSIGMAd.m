%--------------------------------------------------------------------------
% File Name: stylized_rouwenhorst_cSIGMAd.m
% Author: Philip Coyle
% Date Created: 01/08/2019
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/AR1/Rouwenhorst
% stylized_rouwenhorst_cSIGMAd
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
cRHO = 0.8; 
% cSIGMAd_grid = [0.2/100, 0.355/100];
cSIGMAd_grid = [0.15/100, 0.3/100, 0.355/100];

% colors = {'-.b','-k'};
% colors_line = {'b','k'};
colors = {':r','-.b','-k'};
colors_line = {':r','b','k'};

for k = length(cSIGMAd_grid):-1:1
    cSIGMAd = cSIGMAd_grid(k);
    
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd];

    numpts = 31;

    %% Define Rouwenhorst Discretizing State-Space
    [R.del_grid,R.trans_mat,R.s]=rouwenhorst(cRHO,cSIGMAd,numpts);


    %% Target Regime
    % Get initial matrix and solution vector
    [A_tr, b_tr, out_tr] = eqmmat_tr(params,R,numpts);

    y_tr = out_tr(1:numpts);
    pi_tr = out_tr(numpts+1:2*numpts);
    i_tr = out_tr(2*numpts+1:3*numpts);


    % Refine Matrix to account for ZLB
    if sum(i_tr >= 0) == numpts
        converged = 1;
        X_tr = [100*y_tr,400*pi_tr,400*i_tr];
        zlb = find(i_tr == 0,1);
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

            zlb = find(i_tr == 0,1);        
            X_tr = [100*y_tr,400*pi_tr,400*i_tr];
        end

        A_tr = A_up;
        b_tr = b_up;
    end

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
        zlb = find(i_dr == 0,1);
        X_dr = [100*y_dr,400*pi_dr,400*i_dr];
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

            zlb = find(i_dr == 0,1);
            X_dr = [100*y_dr,400*pi_dr,400*i_dr];
        end

        A_dr = A_up;
        b_dr = b_up;
    end

    %% Plotting
    header_tr = {'Output (TR)','Inflation (TR)','Nominal Interest Rate (TR)'};
    header_dr = {'Output (DR)','Inflation (DR)','Nominal Interest Rate (DR)'};
    tr_ylim = [-5 4; -2 1; 0 4];
    dr_ylim = [-6 4; -4 1; 0 3];

    fig(1) = figure(1);
    for i = 1:length(header_tr)
        subplot(2,3,i)
        box on
        hold on
        grid on
        plot(R.del_grid,X_tr(:,i),colors{k},'LineWidth',2)        
        if ~isempty(zlb) && k == length(cSIGMAd_grid)
            line([R.del_grid(zlb) R.del_grid(zlb)], get(gca,'Ylim'),'Color',colors_line{k},'LineWidth',0.5)     
        end
        xlabel('r^n','FontSize',15)
        title(header_tr{i},'FontSize',15,'FontWeight','Normal')
        set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'Ylim',tr_ylim(i,:),'FontSize',15)

        subplot(2,3,i+3)
        box on
        hold on
        grid on
        plot(R.del_grid,X_dr(:,i),colors{k},'LineWidth',2)
        if ~isempty(zlb) && k == length(cSIGMAd_grid)
            line([R.del_grid(zlb) R.del_grid(zlb)], get(gca,'Ylim'),'Color',colors_line{k},'LineWidth',0.5)     
        end
        xlabel('r^n','FontSize',15)
        title(header_dr{i},'FontSize',15,'FontWeight','Normal')
        set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'Ylim',dr_ylim(i,:),'FontSize',15)
    end
end
set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','rouwenhorst_pfs_cSIGMAd.eps');

