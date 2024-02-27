%--------------------------------------------------------------------------
% filename     : CleanData.m
% author       : Philip Coyle
% date created : 02/19/2024
%--------------------------------------------------------------------------
clear all
close all
clc

%% Housekeeping
datapath = '../../Data/';
head = {'Policy Rate', 'Headline Inf.', 'Core Inf.', 'BOJ Core Inf.', 'Core-Core Inf.', 'GDP Defl.'};

savedir = cd;
if ispc 
    savedir = fullfile(savedir, '\..');
    savedir = strcat(savedir,'Model\Final\');
else
    savedir = fullfile(savedir, '/..');
    savedir = strcat(savedir,'/Model/Final/');
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

%% Summary Statistics
autocorr = diag(corr(data(2:end,2:end), data(1:end-1,2:end)));


inflation.mu = round(mean([inf_headline inf_core inf_boj_core inf_core_core defl]),3);
inflation.p50 = round(median([inf_headline inf_core inf_boj_core inf_core_core defl]),3);
inflation.sig2 = round(var([inf_headline inf_core inf_boj_core inf_core_core defl]),3);
inflation.sig = round(std([inf_headline inf_core inf_boj_core inf_core_core defl]),3);
inflation.fracpos = 100*round(mean([inf_headline inf_core inf_boj_core inf_core_core defl] > 0),3);

InflationTable = table(inflation.mu', ...
                       inflation.p50', ...
                       inflation.sig', ...
                       inflation.fracpos', ...
                       'VariableNames', ...
                       ["Mean",...
                        "Median",...
                        "St. Dev",...
                        "Positve (%)"],...
                        'RowNames',...
                        ["Headline Inflation",...
                        "Core Inflation",...
                        "BOJ-Core Inflation",...
                        "Core-Core Inflation",...
                        "GDP Deflator"]);

inx = find(time == 2010);
recent_inflation.mu = round(mean([inf_headline(inx:end) inf_core(inx:end) inf_boj_core(inx:end) inf_core_core(inx:end) defl(inx:end)]),3);
recent_inflation.p50 = round(median([inf_headline(inx:end) inf_core(inx:end) inf_boj_core(inx:end) inf_core_core(inx:end) defl(inx:end)]),3);
recent_inflation.sig2 = round(var([inf_headline(inx:end) inf_core(inx:end) inf_boj_core(inx:end) inf_core_core(inx:end) defl(inx:end)]),3);
recent_inflation.sig = round(std([inf_headline(inx:end) inf_core(inx:end) inf_boj_core(inx:end) inf_core_core(inx:end) defl(inx:end)]),3);
recent_inflation.fracpos = 100*round(mean([inf_headline(inx:end) inf_core(inx:end) inf_boj_core(inx:end) inf_core_core(inx:end) defl(inx:end)] > 0),3);

InflationRecentTable = table(recent_inflation.mu', ...
                       recent_inflation.p50', ...
                       recent_inflation.sig', ...
                       recent_inflation.fracpos', ...
                       'VariableNames', ...
                       ["Mean",...
                        "Median",...
                        "St. Dev",...
                        "Positve (%)"],...
                        'RowNames',...
                        ["Headline Inflation",...
                        "Core Inflation",...
                        "BOJ-Core Inflation",...
                        "Core-Core Inflation",...
                        "GDP Deflator"]);

outputgap.mu = round(mean([gap_boj gap_cabinet]), 3);
outputgap.p50 = round(median([gap_boj gap_cabinet]), 3);
outputgap.sig2 = round(var([gap_boj gap_cabinet]), 3);
outputgap.sig = round(std([gap_boj gap_cabinet]), 3);

OutputGapTable = table(outputgap.mu', ...
                       outputgap.p50', ...
                       outputgap.sig', ...
                       'VariableNames', ...
                       ["Mean",...
                        "Median",...
                        "St. Dev"],...
                        'RowNames',...
                        ["BOJ Measure",...
                        "Cabinet Measure"]);

polrate.mu = round(mean(int),3);
polrate.p50 = round(median(int),3);
polrate.sig2 = round(var(int),3);
polrate.sig = round(std(int),3);

PolRateTable = table(polrate.mu', ...
                     polrate.p50', ...
                     polrate.sig', ...
                     'VariableNames', ...
                     ["Mean",...
                      "Median",...
                      "St. Dev"]);

maketable(InflationTable,"InflationSumStat");
maketable(InflationRecentTable,"InflationRecentSumStat");
maketable(OutputGapTable,"OutputGapSumStat");
maketable(PolRateTable,"PolRateSumStat");

%% Time Series Plots
XX = [int, inf_headline, inf_core, inf_boj_core, inf_core_core, defl];

fig(1) = figure(1);
for i = 1:length(head)
    subplot(3,2,i)
    box on
    grid on
    hold on
    plot(time, XX(:,i), 'LineWidth', 2, 'color', 'k')
    

    set(gca,'fontsize',10,'XLim',[1998 2019])
    yline(0,'k','linewidth',1)
    xlabel('Time','fontsize',15)
    title(head{i},'fontsize',15, 'FontWeight','normal')
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-dpdf',strcat( savedir,'TSplots.pdf'));

%% Rolling Average Plots (Rolling End)
endinx = find(time == 2008);
XX_rolling_end = zeros(length(time) - endinx + 1, 6);
for i = 1:size(XX_rolling_end,1)
    XX_rolling_end(i,:) = mean(XX(1:endinx + i - 1, :));
end

fig(2) = figure(2);
for i = 1:length(head)
    subplot(3,2,i)
    box on
    grid on
    hold on
    plot(time(endinx:end), XX_rolling_end(:,i), 'LineWidth', 2, 'color', 'k')
    

    set(gca,'fontsize',10,'XLim',[2008 2019])
    yline(0,'k','linewidth',1)
    xlabel('Time','fontsize',15)
    title(head{i},'fontsize',15, 'FontWeight','normal')
end

set(fig(2),'PaperOrientation','Landscape');
set(fig(2),'PaperPosition',[0 0 11 8.5]);
print(fig(2),'-dpdf',strcat( savedir,'RollingMeanEnd.pdf'));


%% Rolling Average Plots (Rolling Start)
startinx = find(time == 2013);
XX_rolling_start = zeros(length(time) - startinx + 1, 6);
for i = 1:startinx
    XX_rolling_start(i,:) = mean(XX(i:end, :));
end

fig(3) = figure(3);
for i = 1:length(head)
    subplot(3,2,i)
    box on
    grid on
    hold on
    plot(time(1:startinx), XX_rolling_start(:,i), 'LineWidth', 2, 'color', 'k')
    

    set(gca,'fontsize',10,'XLim',[1998 2013])
    yline(0,'k','linewidth',1)
    xlabel('Time','fontsize',15)
    title(head{i},'fontsize',15, 'FontWeight','normal')
end

set(fig(3),'PaperOrientation','Landscape');
set(fig(3),'PaperPosition',[0 0 11 8.5]);
print(fig(3),'-dpdf',strcat( savedir,'RollingMeanStart.pdf'));


%% Kernel Density 
for j = 1:3
    if j == 1
        startingyear = 1998;
    elseif j == 2
        startingyear = 2005;
    elseif j == 3
        startingyear = 2010;
    end
    startinx = find(time == startingyear);
    
    for i = 1:6
        [f, xi] = ksdensity(XX(startinx:end,i));
        XX_kdensity(:,i) = f;
        xx_kdensity(:,i) = xi;
    end

    fig(3+j) = figure(3+j);
    for i = 1:length(head)
        subplot(3,2,i)
        box on
        grid on
        hold on
        plot(xx_kdensity(:,i), XX_kdensity(:,i), 'LineWidth', 2, 'color', 'k')
        
    
        set(gca,'fontsize',10)
        ylabel('Density','fontsize',15)
        title(head{i},'fontsize',15, 'FontWeight','normal')
    end
    
    
    set(fig(3+j),'PaperOrientation','Landscape');
    set(fig(3+j),'PaperPosition',[0 0 11 8.5]);
    print(fig(3+j),'-dpdf',strcat( savedir,'KernDensities_',num2str(startingyear),'.pdf'));
end