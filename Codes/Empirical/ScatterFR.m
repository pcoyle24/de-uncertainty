%--------------------------------------------------------------------------
% File Name: ScatterFR.m
% Author: Philip Coyle
% Date Created: 02/26/2024
%--------------------------------------------------------------------------

clear all
close all
clc

%% Housekeeping
datapath = '../../Data/';
head = {'Headline Inflation', 'Core Inflation', 'BOJ Core Inflation', 'Core-Core Inflation', 'GDP Deflator'};

savedir = cd;
if ispc 
    savedir = fullfile(savedir, '\..');
    savedir = strcat(savedir,'Model\Final\');
else
    savedir = fullfile(savedir, '/..');
    savedir = strcat(savedir,'/Model/Final/');
end


%% Paramaters
cSIGMA = 1.05;
cKAPPA = 0.04;
cPHIpi = 1.75;
cRstar = 0.5/400;
cBET = 1/(1 + cRstar);
cRHO = 0.81; 
cPItarg1= 1/400;
cPItarg2= 2/400;
cSIGMAd = 0.15/100;

params = [cBET;cSIGMA;cKAPPA;cPHIpi;cRstar;cSIGMAd];

nstate = 21;
% Find Middle State
mid = (nstate+1)/2;

% Set Grid Intervals
pi_m_low = -2.5/400;
pi_m_high = cPItarg2 + 0.5/400;
pts = 1001;
pi_m = linspace(pi_m_low,pi_m_high, pts)';

% Taylor Rule
tr1 = max(0,cRstar + cPItarg1 + cPHIpi*(pi_m - cPItarg1));
tr2 = max(0,cRstar + cPItarg2 + cPHIpi*(pi_m - cPItarg2));

% 'Riskless' Fisher Relation
fr = pi_m + cRstar;

%% Get DSS
[~,inx2] = sort(abs(tr2 - fr));
DSS2 = inx2(1:2);
it = 2;
while abs(DSS2(1) - DSS2(2)) <= 10
    DSS2(2) = inx2(it + 1);
    it = it + 1;
end

[~,inx1] = sort(abs(tr1 - fr));
DSS1 = inx1(1:2);
it = 2;
while abs(DSS1(1) - DSS1(2)) <= 10
    DSS1(2) = inx1(it + 1);
    it = it + 1;
end

%% Read in data
% values = {'A2:I104'}; % from 1998Q1 to 2023Q3
values = {'A2:I89'}; % from 1998Q1 to 2019Q2
file_name = strcat(datapath,'cleaned/data.xlsx');
data = xlsread(file_name,char(values));

time = data(:,1);
quarter = 4*mod(time,1) + 1;
inf_headline = data(:,2);
inf_core = data(:,3);
inf_boj_core = data(:,4);
inf_core_core = data(:,5);
gap_boj = data(:,6);
gap_cabinet = data(:,7);
defl = data(:,8);
int = data(:,9);


startinx = find(time == 2013);
XX = [inf_headline, inf_core, inf_boj_core, inf_core_core, defl];
%% Plotting
for i = 1:size(XX,2)
    fig(i) = figure(i);
    grid on
    box on
    hold on
    
    % Plot Lines
    h(1) = plot(400*pi_m, 400*tr2,'Color','k','LineStyle','--','LineWidth',2);
    h(2) = plot(400*pi_m, 400*tr1,'Color','k','LineStyle','-.','LineWidth',2);
    h(3) = plot(400*pi_m, 400*fr,'Color','k','LineStyle','-','LineWidth',2);
    
%     % Plot DSS
%     pdss1 = plot(400*pi_m(DSS(1)),400*tr(DSS(1)));
%     pdss2 = plot(400*pi_m(DSS(2)),400*tr(DSS(2))); 
%     set(pdss1,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
%     set(pdss2,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
    
    % Plot Data
    d(1) = scatter(XX(1:startinx - 1,i),int(1:startinx - 1),'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
    d(2) = scatter(XX(startinx:end,i),int(startinx:end),'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
    
    
    % Clean up x and y axis
    set(gca,'XLim',[400*pi_m(1), 400*pi_m(end)],'YLim',[-1, 4],'FontSize',25)
    xlabel('Inflation Rate')
    ylabel('Nominal Interest Rate')
    title(head{i}, 'FontWeight','normal')
    
    % Make legend
    L = legend([h(1) h(2) h(3) d(1) d(2)],'Taylor Rule (\pi* = 2%)','Taylor Rule (\pi* = 1%)','Standard Fisher Relation', 'Data: 1998 - 2012','Data: 2013 - 2019');
    set(L,'Location','NorthWest','Fontsize',20)
    
    set(fig(i),'PaperOrientation','Landscape');
    set(fig(i),'PaperPosition',[0 0 11 8.5]);
    if i == 1
        print(fig(i),'-dpdf',strcat(savedir,'scatter_headlineinf_FR.pdf'));
    elseif i == 2
        print(fig(i),'-dpdf',strcat(savedir,'scatter_coreinf_FR.pdf'));
    elseif i == 3
        print(fig(i),'-dpdf',strcat(savedir,'scatter_bojcoreinf_FR.pdf'));
    elseif i == 4
        print(fig(i),'-dpdf',strcat(savedir,'scatter_corecoreinf_FR.pdf'));
    elseif i == 5
        print(fig(i),'-dpdf',strcat(savedir,'scatter_gdpdefl_FR.pdf'));
    end
end
