%--------------------------------------------------------------------------
% File Name: stylized_rouwenhorst_tr.m
% Author: Philip Coyle
% Date Created: 01/08/2019
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/AR1/Rouwenhorst
% stylized_rouwenhorst_tr
%--------------------------------------------------------------------------

clear all
close all
clc

warning('off')

%% Paramaters
cBET = 1/1.0025;
cSIGMA = 1;
cKAPPA = 0.02;
cPHIpi = 2;
cRstar = 1/400;
cRHO = 0.8; 
cSIGMAd = 0.2/100;

params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cRHO;cSIGMAd];

numpts = 11;

%% Define Rouwenhorst Discretizing State-Space
[R.del_grid,R.trans_mat,R.s]=rouwenhorst(cRHO,cSIGMAd,numpts);

%% Main Code
% Get initial matrix and solution vector
[A, b, out] = eqmmat_tr(params,R,numpts);

y_out = out(1:numpts);
pi_out = out(numpts+1:2*numpts);
i_out = out(2*numpts+1:3*numpts);


% Refine Matrix to account for ZLB
if sum(i_out >= 0) == numpts
    converged = 1;
    X = [100*y_out,400*pi_out,400*i_out];
    zlb = find(i_out == 0,1);
else
    converged = 0;
end

while converged == 0
    [A_up, b_up, out_up] = eqmrefine_tr(params,R,numpts,A,b);
    
    i_out = out_up(2*numpts+1:3*numpts);

    if sum(i_out >= 0) == numpts
        converged =1;
        
        y_out = out_up(1:numpts);
        pi_out = out_up(numpts+1:2*numpts);
        i_out = out_up(2*numpts+1:3*numpts);
        
        % Check equilibirum existence
        i_real = cRstar + cPHIpi*pi_out;
        inx = find(i_real>0);
        if sum(abs(i_out(inx) - i_real(inx)) < eps) ~= length(inx)
            error('The equilibrium does not exist. Reduce shock size.')
        end
        
        zlb = find(i_out == 0,1);        
        X = [100*y_out,400*pi_out,400*i_out];
    end
    
    A = A_up;
    b = b_up;
end

%% Plotting
header = {'Output','Inflation','Nominal Interest Rate'};

fig(1) = figure(1);
for i = 1:length(header)
    subplot(2,3,i)
    box on
    hold on
    grid on
    plot(R.del_grid,X(:,i),'Color','k','LineWidth',2)
%     if ~isempty(zlb)
%         line([R.del_grid(zlb) R.del_grid(zlb)], get(gca,'Ylim'),'Color','b','LineWidth',1) 
%     end
    xlabel('\delta','FontSize',15)
    title(header{i},'FontSize',15,'FontWeight','Normal')
    set(gca,'Xlim',[R.del_grid(1) R.del_grid(end)],'FontSize',15)
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','rouwenhorst_pfs_tr.eps');
