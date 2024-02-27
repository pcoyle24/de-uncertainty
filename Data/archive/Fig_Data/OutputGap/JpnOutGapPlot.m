%% About File
%--------------------------------------------------------------------------
% filename     : JpnOutGapPlot.m
% author       : Philip Coyle
% date created : 08/27/2018
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Data/UpdatedData/OutputGap
% JpnOutGapPlot.m
%--------------------------------------------------------------------------

%% Main Code
clear all
close all
clc

var_names = {'BoJ Output Gap','Cabinet Office Output Gap'};
sheet_names = {'Sheet1'};
period_vec = 1983:0.25:2018;

values = {'B14:C154'}; % from 1983 to 2018
value_mat = xlsread('gap.xlsx',char(sheet_names(1)),char(values));


fig1 = figure(1);
i_sub_mark = 1;
for j = 1:length(var_names)
    subplot(2,2,i_sub_mark)
    box on
    hold on
    grid on
    
    h1 = plot(period_vec',value_mat(:,i_sub_mark),'k','LineWidth',2);
    title(var_names(j),'FontSize',14,'FontWeight','normal')
    set(gca,'FontSize',12,'XLim',[1983 2018],'XTick',1983:7:2018)
    xlabel('Year','FontSize',12)

    set(gca,'YLim',[-8 6],'YTick',-8:2:6)

    i_sub_mark = i_sub_mark + 1;
end

set(fig1,'PaperOrientation','Landscape');
set(fig1,'PaperPosition',[0 0 11 8.5]);
print(fig1,'-depsc','JpnOutGapPlot.eps');

