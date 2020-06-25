%% About File
%--------------------------------------------------------------------------
% filename     : FisherRelation.m
% author       : Philip Coyle
% date created : 08/28/2018
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Codes/StylizedModel/SS_Opt_Inf/FisherRelation
% FisherRelation.m
%--------------------------------------------------------------------------

%% Main Code

clear all
close all
clc

x = 0:1:50;
y1 = 0.5.*x - 3;
y2 = 0.6.*x - 1.5;
y3(1,(1:length(x))) = 0;
y4(1,(1:length(x))) = 0;

for i = 1:length(x)
    if i <= 16
        y3(1,i) = 0;
        y4(1,i) = 0;
    elseif i <= 22
        y3(1,i) = x(1,i) - 16;
        y4(1,i) = 0;
    else
        y3(1,i) = x(1,i) - 16;
        y4(1,i) = x(1,i)*15.125/14.25 - 22*15.125/14.25;       
    end		
end

fig1 = figure(1);
box off
hold on
grid off
h1 = plot(x,y2,'k',x,y3,'r','LineWidth',2);
set(gca,'XLim',[0 50],'XTick',[36.25],'XTickLabel',{'2'},'YLim',[-8 25],'YTick',[0],'YTickLabel',{'1'},'FontSize',22)
plot([36.25 36.25],[20.25 20.25],'.k','MarkerSize',25);
plot([2.5 2.5],[0 0],'.k','MarkerSize',25);
tr = text(34.2,21.1,'TR');
dr = text(3.5,-1.25,'DR');
inf = text(52.5,-8,'\bf \Pi','Interpreter','tex');
R = text(-2,24,'\bf R','Interpreter','tex');
set(inf,'FontWeight','bold','HorizontalAlignment','center','FontSize',30)
set(R,'FontWeight','bold','HorizontalAlignment','center','FontSize',30)
set(tr,'FontWeight','bold','HorizontalAlignment','center','FontSize',22)
set(dr,'FontWeight','bold','HorizontalAlignment','center','FontSize',22)
c = [0.41,0.46];
d = [0.65,0.61];
arrow_2 = annotation('textarrow',c,d,'String',({'Fisher Relation'}));
set(arrow_2,'FontWeight','bold','HorizontalAlignment','center','FontSize',22)
e = [0.43,0.41];
f = [0.26,0.34];
arrow_3 = annotation('textarrow',e,f,'String',{'Taylor Rule'});
set(arrow_3,'Color','r','FontWeight','bold','HorizontalAlignment','center','FontSize',22)
annotation('arrow',[0.13 .92],[0.11 0.11])
annotation('arrow',[0.13 0.13],[0.11 0.94])


set(fig1,'PaperOrientation','Landscape');
set(fig1,'PaperPosition',[0 0 11 8.5]);
print(fig1,'-depsc','FisherRelation.eps');
