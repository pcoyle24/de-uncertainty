%--------------------------------------------------------------------------
% File Name: PFs_cPHIpi_h.m
% Author: Philip Coyle
% Date Created: 01/08/2019
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
cPItarg= 2/400; 
cSIGMAd_grid = [ 0.55/100,  0.60/100,  cSIGMAd_max.h-eps];

colors = {':r','-.b','-k'};
colors_line = {':r','b','k'};

for k = length(cSIGMAd_grid):-1:1
    cSIGMAd = cSIGMAd_grid(k);
    
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd;cPItarg];

    numpts = 21;

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
        zlb_tr = find(i_tr < 1e-10,1);
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
            i_real_tr = cRstar + cPItarg + cPHIpi*(pi_tr - cPItarg);
            inx = find(i_real_tr>0);
            if sum(abs(i_tr(inx) - i_real_tr(inx)) < eps) ~= length(inx)
                error('The equilibrium does not exist. Reduce shock size.')
            end

            zlb_tr = find(i_tr < 1e-10,1);  
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

    i_real_dr = cRstar + cPItarg + cPHIpi*(pi_dr - cPItarg);

    % Refine Matrix to account for ZLB
    if sum(i_dr > 0) == sum(i_real_dr > 0)
        converged = 1;
        zlb_dr = find(i_dr < 1e-10,1);
        X_dr = [100*y_dr,400*pi_dr,400*i_dr];
    else
        converged = 0;
    end
    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_dr(params,R,numpts,A_dr,b_dr);

        pi_dr = out_up(numpts+1:2*numpts);
        i_dr = out_up(2*numpts+1:3*numpts);
        i_real_dr = cRstar + cPItarg + cPHIpi*(pi_dr - cPItarg);


        if sum(i_dr > 0) == sum(i_real_dr > 0)
            converged =1;

            y_dr = out_up(1:numpts);
            pi_dr = out_up(numpts+1:2*numpts);
            i_dr = out_up(2*numpts+1:3*numpts);

            % Check equilibirum existence
            if sum(i_dr<0) > 0
                error('The equilibrium does not exist. Reduce shock size.');
            end

            zlb_dr = find(i_dr < 1e-10,1);
            X_dr = [100*y_dr,400*pi_dr,400*i_dr];
        end

        A_dr = A_up;
        b_dr = b_up;
    end

    %% Plotting
    header_tr = {'Output (TE)','Inflation (TE)','Policy Rate (TE)'};
    header_dr = {'Output (DE)','Inflation (DE)','Policy Rate (DE)'};
    tr_ylim = [-30 10; -10 6; 0 20];
    dr_ylim = [-30 10; -10 6; 0 20];

    fig(1) = figure(1);
    for i = 1:length(header_tr)
        subplot(2,3,i)
        box on
        hold on
        grid on
        plot(R.del_grid,X_tr(:,i),colors{k},'LineWidth',2)        
        if ~isempty(zlb_tr) && k == length(cSIGMAd_grid)
            line([R.del_grid(zlb_tr) R.del_grid(zlb_tr)], get(gca,'Ylim'),'Color',colors_line{k},'LineWidth',0.5)     
        end
        xlabel('\delta','FontSize',20) 
        title(header_tr{i},'FontSize',20,'FontWeight','Normal')
        set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'Ylim',tr_ylim(i,:),'FontSize',20)

        subplot(2,3,i+3)
        box on
        hold on
        grid on
        if i == 3
            h(k) = plot(R.del_grid,X_dr(:,i),colors{k},'LineWidth',2);
        else
           plot(R.del_grid,X_dr(:,i),colors{k},'LineWidth',2);
        end
        if ~isempty(zlb_dr) && k == length(cSIGMAd_grid)
            line([R.del_grid(zlb_dr) R.del_grid(zlb_dr)], get(gca,'Ylim'),'Color',colors_line{k},'LineWidth',0.5)     
        end
        xlabel('\delta','FontSize',20) 
        title(header_dr{i},'FontSize',20,'FontWeight','Normal')
        set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'Ylim',dr_ylim(i,:),'FontSize',20)
    end
end
L = legend([h(1) h(2) h(3)],'\sigma_{\epsilon} = 0.55/100','\sigma_{\epsilon} = 0.60/100','\sigma_{\epsilon} = \sigma_{\epsilon}^{max}');
set(L,'Location','NorthEast','Fontsize',14)
% set(L,'Orientation','horizontal','Location','SouthOutside','Fontsize',18)
% lp = get(L,'position');  
% set(L,'position',[-0.375*lp(1) lp(2)/1.6 lp(3) lp(4)]);  % sets the legend farther down
% legend('boxoff')

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
print(fig(1),'-dpdf',strcat(savedir,'PFs_cPHIpi_h.pdf'));


