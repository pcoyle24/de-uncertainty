%--------------------------------------------------------------------------
% File Name: stylized_rouwenhorst_moments_cPHIpi_m.m
% Author: Philip Coyle
% Date Created: 01/18/2019
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Uncertainty/Draft/Figs/ar1/nstate/moments
% matlab -nodesktop -nosplash -r stylized_rouwenhorst_moments_cPHIpi_m
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
cSIGMAd_grid = linspace(0,0.4/100,101);

numpts = 21;

%% Housekeeping
pf_y_dr = nan(numpts, length(cSIGMAd_grid));
pf_pi_dr = nan(numpts, length(cSIGMAd_grid));
pf_i_dr = nan(numpts, length(cSIGMAd_grid));

mu_y_dr = nan(length(cSIGMAd_grid),1);
mu_ym_dr = nan(length(cSIGMAd_grid),1);
sig_y_dr = nan(length(cSIGMAd_grid),1);
sig_ym_dr = nan(length(cSIGMAd_grid),1);

mu_pi_dr = nan(length(cSIGMAd_grid),1);
mu_pim_dr = nan(length(cSIGMAd_grid),1);
sig_pi_dr = nan(length(cSIGMAd_grid),1);
sig_pim_dr = nan(length(cSIGMAd_grid),1);

mu_i_dr = nan(length(cSIGMAd_grid),1);
mu_im_dr = nan(length(cSIGMAd_grid),1);
sig_i_dr = nan(length(cSIGMAd_grid),1);
sig_im_dr = nan(length(cSIGMAd_grid),1);

mu_elb_dr = nan(length(cSIGMAd_grid),1);
mu_elbm_dr = nan(length(cSIGMAd_grid),1);
sig_elb_dr = nan(length(cSIGMAd_grid),1);
sig_elbm_dr = nan(length(cSIGMAd_grid),1);


tic;
for i = 1:length(cSIGMAd_grid)
    cSIGMAd = cSIGMAd_grid(i);
    params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd;cPItarg];


    %% Define Rouwenhorst Discretizing State-Space
    [R.del_grid,R.trans_mat,R.s]=rouwenhorst(cRHO,cSIGMAd,numpts);
    
    if i == 1
        rss_inx = (numpts+1)/2;
    else
        rss_inx = find(abs(R.del_grid - 0) < 1e-10);
    end

    R.m = R.trans_mat(rss_inx,:)';

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
    else
        converged = 0;
    end
    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_dr(params,R,numpts,A_dr,b_dr);

        pi_dr = out_up(numpts+1:2*numpts);
        i_dr = out_up(2*numpts+1:3*numpts);
        i_real_dr = cRstar + cPItarg + cPHIpi*(pi_dr - cPItarg);


        if sum(i_dr > 0) == sum(i_real_dr > 0)
            converged = 1;

            y_dr = out_up(1:numpts);
            pi_dr = out_up(numpts+1:2*numpts);
            i_dr = out_up(2*numpts+1:3*numpts);

            % Check equilibirum existence
            if sum(i_dr<0) > 0
%                 disp('The equilibrium does not exist. Reduce shock size.')
                converged = -1;
            end
        end

        A_dr = A_up;
        b_dr = b_up;
    end
    
    if converged == 1

        % Store Policy Functions
        pf_y_dr(:,i) = 100*y_dr;
        pf_pi_dr(:,i) =400*pi_dr;
        pf_i_dr(:,i) = 400*i_dr;

        % Calculate Moments with stationary distribution R.s
        % Get Expected Values
        mu_y_dr(i) = 100*sum(R.s.*y_dr);
        mu_pi_dr(i) = 400*sum(R.s.*pi_dr);
        mu_i_dr(i) = 400*sum(R.s.*i_dr);
        mu_elb_dr(i) = 100*sum(R.s.*(i_dr == 0));

        % Get Standard Deviations
        sig_y_dr(i) = 100*real(sqrt(sum(R.s.*(y_dr).^2) - sum(R.s.*y_dr)^2));
        sig_pi_dr(i) = 400*real(sqrt(sum(R.s.*(pi_dr).^2) - sum(R.s.*pi_dr)^2));
        sig_i_dr(i) = 400*real(sqrt(sum(R.s.*(i_dr).^2) - sum(R.s.*i_dr)^2));
        sig_elb_dr(i) = 100*real(sqrt(sum(R.s.*((i_dr == 0)).^2) - sum(R.s.*(i_dr == 0))^2));


%         % Calculate Moments with stationary distribution R.s
%         % Get Expected Values
%         mu_ym_dr(i) = sum(R.m.*y_dr*100);
%         mu_pim_dr(i) = sum(R.m.*pi_dr*400);
%         mu_im_dr(i) = sum(R.m.*i_dr*400);
%         mu_elbm_dr(i) = sum(R.m.*(i_dr == 0)*100);
% 
%         % Get Standard Deviations
%         sig_ym_dr(i) = real(sqrt(sum(R.m.*(y_dr*100).^2) - sum(R.m.*y_dr*100)^2));
%         sig_pim_dr(i) = real(sqrt(sum(R.m.*(pi_dr*400).^2) - sum(R.m.*pi_dr*400)^2));
%         sig_im_dr(i) = real(sqrt(sum(R.m.*(i_dr*400).^2) - sum(R.m.*i_dr*400)^2));
%         sig_elbm_dr(i) = real(sqrt(sum(R.m.*((i_dr == 0)*100).^2) - sum(R.m.*(i_dr == 0)*100)^2));

    end
    
    if mod(i,10) == 0
        disp(strcat('cSIGMAd = '," ",num2str(i)," ",'of'," ",num2str(length(cSIGMAd_grid))))
    end
end
toc;

%% Plotting
X_dr = [pf_y_dr(rss_inx,:)', pf_pi_dr(rss_inx,:)', pf_i_dr(rss_inx,:)',...
        mu_y_dr, mu_pi_dr, mu_i_dr,...
        sig_y_dr, sig_pi_dr, sig_i_dr];

header = {'y_{RSS}','\pi_{RSS}','i_{RSS}',...
          'E[y]','E[\pi]','E[i]',...
          'StDev[y]','StDev[\pi]','StDev[i]'};

fig(1) = figure(1);
for i = 1:length(header)
    subplot(3,3,i)
    box on
    grid on
    hold on
    if i == length(header)
        h(1) = plot(cSIGMAd_grid,X_dr(:,i),'Color','k','LineStyle','-','LineWidth',2);
    else
        plot(cSIGMAd_grid,X_dr(:,i),'Color','k','LineStyle','-','LineWidth',2)
    end
    
    xline(0.125/100, 'LineWidth',1,'Color','r')
    if i == 4
        yline(-1.11, 'LineWidth',1,'Color','b')
        yline(-0.57, 'LineWidth',1,'Color','b')
    elseif i == 5
        yline(0.01, 'LineWidth',1,'Color','b')
        yline(0.37, 'LineWidth',1,'Color','b')
    elseif i == 6
        yline(0.0, 'LineWidth',1,'Color','b')
        yline(0.02, 'LineWidth',1,'Color','b')
    elseif i == 7
        yline(1.49, 'LineWidth',1,'Color','b')
        yline(1.87, 'LineWidth',1,'Color','b')
    elseif i == 8
        yline(1.41, 'LineWidth',1,'Color','b')
        yline(1.78, 'LineWidth',1,'Color','b')
    elseif i == 9
        yline(0.1, 'LineWidth',1,'Color','b')
        yline(0.2, 'LineWidth',1,'Color','b')
    end
    xlabel('\sigma_{\epsilon}','FontSize',15)
    title(header{i},'FontSize',15,'FontWeight','Normal')
    set(gca,'Xlim',[cSIGMAd_grid(1) cSIGMAd_grid(end)],'FontSize',15)
end


% pltinx = find(isnan(mu_y_dr),1) - 1;
% if isempty(pltinx)
%     pltinx = length(R.del_grid);
% end
% figure(2);
% subplot(2,3,1)
% box on
% grid on
% hold on
% plot(R.del_grid, pf_y_dr(:,pltinx), 'LineWidth',2)
% xline(0, 'LineWidth',1)
% set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'FontSize',15)
% xlabel('\delta','FontSize',15)
% title('Output Gap','FontSize',15,'FontWeight','Normal')
% 
% subplot(2,3,2)
% box on
% grid on
% hold on
% plot(R.del_grid, pf_pi_dr(:,pltinx), 'LineWidth',2)
% xline(0, 'LineWidth',1)
% set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'FontSize',15)
% xlabel('\delta','FontSize',15)
% title('Inflation','FontSize',15,'FontWeight','Normal')
% 
% subplot(2,3,3)
% box on
% grid on
% hold on
% plot(R.del_grid, pf_i_dr(:,pltinx), 'LineWidth',2)
% xline(0, 'LineWidth',1)
% set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'FontSize',15)
% xlabel('\delta','FontSize',15)
% title('Policy Rate','FontSize',15,'FontWeight','Normal')

savedir = cd;
if ispc 
    savedir = fullfile(savedir, '..\..\..');
    savedir = strcat(savedir,'\Final\');
else
    savedir = fullfile(savedir, '../../..');
    savedir = strcat(savedir,'/Final/');
end

% set(fig(1),'PaperOrientation','Landscape');
% set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-dpdf',strcat(savedir,'moment_matching.pdf'));