%% About File
%--------------------------------------------------------------------------
% filename     : JpnTimeSeries.m
% author       : Philip Coyle
% date created : 08/27/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_1
% JpnTimeSeries.m
%--------------------------------------------------------------------------

%% Main Code
clear all
close all
clc

sheet_names = {'JpnData'};
var_names = {'Inflation (Annualized %)','Output Gap (%)','Policy Rate (Annualized %)'};
values = {'B66:D194'}; % from 1983Q1 to 2017Q1
period_vec = 1986:0.25:2018;
value_mat = xlsread('jpndata.xlsx',char(sheet_names(1)),char(values));


fig1 = figure(1);
i_sub_mark = 1;
for j = 1:length(var_names)
    subplot(2,3,i_sub_mark)
    box on
    hold on
    grid on
%     if i_sub_mark == 3
%         A1 = area([2001 2001.75],[4.95 4.95],-9.95,'LineStyle','none');
%         A2 = area([2007.75 2009.25],[4.95 4.95],-9.95,'LineStyle','none');
%         set([A1 A2],'FaceColor',[0.85 0.85 0.85]) % grey color
%     elseif i_sub_mark == 2
%         A1 = area([2001 2001.75],[2.99 2.99],-0.99,'LineStyle','none');
%         A2 = area([2007.75 2009.25],[2.99 2.99],-0.99,'LineStyle','none');
%         set([A1 A2],'FaceColor',[0.85 0.85 0.85]) % grey color
%     else
%         A1 = area([2001 2001.75],[7.98 7.98],0.02,'LineStyle','none');
%         A2 = area([2007.75 2009.25],[7.98 7.98],0.02,'LineStyle','none');
%         set([A1 A2],'FaceColor',[0.85 0.85 0.85]) % grey color  
%     else 
%         A1 = area([2001 2001.75],[3.99 3.99],1.51,'LineStyle','none');
%         A2 = area([2007.75 2009.25],[3.99 3.99],1.51,'LineStyle','none');
%         set([A1 A2],'FaceColor',[0.85 0.85 0.85]) % grey color  
%     end
    
    h1 = plot(period_vec',value_mat(:,i_sub_mark),'k','LineWidth',2);
%     h2 = line([2009 2009],get(gca,'YLim'),'Color','k','LineStyle','--','LineWidth',1.5);
%     line([2016 2016],get(gca,'YLim'),'Color','k','LineStyle','--','LineWidth',1.5);
    title(var_names(j),'FontSize',14,'FontWeight','normal')
    set(gca,'FontSize',12,'XLim',[1986 2018],'XTick',1986:8:2018)
    xlabel('Year','FontSize',12)
    if i_sub_mark == 3
        set(gca,'YLim',[-1 10],'YTick',0:2:10)
        h3 = line(get(gca,'XLim'),[0 0],'Color','k','LineStyle','-','LineWidth',1);
    elseif i_sub_mark == 2
        set(gca,'YLim',[-7 6],'YTick',-6:3:6)
        h3 = line(get(gca,'XLim'),[0 0],'Color','k','LineStyle','-','LineWidth',1);
    elseif i_sub_mark == 1
        set(gca,'YLim',[-2 4],'YTick',-2:2:4)
%         h2 = line([2009 2009],get(gca,'YLim'),'Color','k','LineStyle','-','LineWidth',1);
        h3 = line(get(gca,'XLim'),[0 0],'Color','k','LineStyle','-','LineWidth',1);
    end
    i_sub_mark = i_sub_mark + 1;
end

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig1,'PaperOrientation','Landscape');
set(fig1,'PaperPosition',[0 0 11 8.5]);
print(fig1,'-depsc','JpnTimeSeries.eps');
print(fig1,'-depsc',strcat(savedir,'JpnTimeSeries.eps'));

