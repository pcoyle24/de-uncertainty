%% About File
%--------------------------------------------------------------------------
% filename     : JpnNomIntPlot.m
% author       : Philip Coyle
% date created : 08/27/2018
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Data/UpdatedData/IntRate
% JpnNomIntPlot.m
%--------------------------------------------------------------------------

%% Main Code
clear all
close all
clc

var_names = {'Policy Rate'};
sheet_names = {'NomRate'};
period_vec = 1985.5:0.25:2018.5;

values = {'B26:B158'}; % from 1985Q3 to 2018Q3
value_mat = xlsread('NomRate.xlsx',char(sheet_names(1)),char(values));

fig1 = figure(1);
for j = 1:length(var_names)
    box on
    hold on
    grid on
    
    h1 = plot(period_vec',value_mat,'k','LineWidth',2);
    title(var_names(j),'FontSize',14,'FontWeight','normal')
    set(gca,'FontSize',12,'XLim',[1985.5 2018.5],'XTick',1986:4:2018)
    xlabel('Year','FontSize',12)

    set(gca,'YLim',[-0.1 10],'YTick',0:2:10)

end

set(fig1,'PaperOrientation','Landscape');
set(fig1,'PaperPosition',[0 0 11 8.5]);
print(fig1,'-depsc','JpnNomIntPlot.eps');

