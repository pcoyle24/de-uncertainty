%% About File
%--------------------------------------------------------------------------
% filename     : JpnInflationPlot.m
% author       : Philip Coyle
% date created : 08/27/2018
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Data/UpdatedData/Inflation
% JpnInflationPlot.m
%--------------------------------------------------------------------------

%% Main Code
clear all
close all
clc

var_names = {'GDP Defl (All)','GDP Defl (Private Consumption)','GDP Defl (Consumption of Households)','Core CPI','BoJ-Core CPI','Core-Core CPI'};
sheet_names = {'Sheet1'};
period_vec = 1994.25:0.25:2018.25;


cd GDPDeflator
values_defl = {'F4:H100'}; % from 1994Q2 to 2018Q2
value_mat_defl = xlsread('GDPDefl.xlsx',char(sheet_names(1)),char(values_defl));

cd ../CPI
values_defl = {'B97:E193'}; % from 1994Q2 to 2018Q2
value_mat_cpi = xlsread('AnnInf_q.xlsx',char(sheet_names(1)),char(values_defl));

cd ..
value_mat = [value_mat_defl,value_mat_cpi];

fig1 = figure(1);
i_sub_mark = 1;
for j = 1:length(var_names)
    subplot(2,3,i_sub_mark)
    box on
    hold on
    grid on
    
    h1 = plot(period_vec',value_mat(:,i_sub_mark),'k','LineWidth',2);
    title(var_names(j),'FontSize',14,'FontWeight','normal')
    set(gca,'FontSize',12,'XLim',[1994 2018],'XTick',1994:6:2018)
    xlabel('Year','FontSize',12)

    set(gca,'YLim',[-5 10],'YTick',-5:5:10)

    i_sub_mark = i_sub_mark + 1;
end

set(fig1,'PaperOrientation','Landscape');
set(fig1,'PaperPosition',[0 0 11 8.5]);
print(fig1,'-depsc','JpnInflationPlot.eps');

