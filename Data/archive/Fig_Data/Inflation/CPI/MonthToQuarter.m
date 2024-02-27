%% About File
%--------------------------------------------------------------------------
% filename     : MonthToQuarter.m
% author       : Philip Coyle
% date created : 08/29/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Codes/StylizedModel/SS_Opt_Inf/Data
% MonthToQuarter.m
%--------------------------------------------------------------------------

clear all
close all
clc

% Housekeeping
sheet_names = {'Sheet1'};
values = {'B2:E584'}; % from 1983Q1 to 2017Q1
cpi_m = xlsread('CPI_m.xlsx',char(sheet_names(1)),char(values));

cpi_headline_m = cpi_m(:,1);
cpi_core_m = cpi_m(:,2);
cpi_boj_core_m = cpi_m(:,3);
cpi_core_core_m = cpi_m(:,4);

quart = 1970:.25:2018.25;
cpi_headline_q = zeros(length(quart),1);
cpi_core_q = zeros(length(quart),1);
cpi_boj_core_q = zeros(length(quart),1);
cpi_core_core_q = zeros(length(quart),1);

% Convert Monthly Series to Quarterly Series
j = 1;
for i = 1:3:length(cpi_m)-1
    cpi_headline_q(j) = sum(cpi_headline_m(i:i+2))/3;
    cpi_core_q(j) = sum(cpi_core_m(i:i+2))/3;
    cpi_boj_core_q(j) = sum(cpi_boj_core_m(i:i+2))/3;
    cpi_core_core_q(j) = sum(cpi_core_core_m(i:i+2))/3;
    
    j = j+1;
end

% Write to Excel
filename = 'CPI_q.xlsx';
cpi_q = [quart',cpi_headline_q,cpi_core_q,cpi_boj_core_q,cpi_core_core_q];
xlswrite(filename,cpi_q)

% Annaulized Log-Difference Calculation
inf_headline_q = zeros(length(quart)-1,1);
inf_core_q = zeros(length(quart)-1,1);
inf_boj_core_q = zeros(length(quart)-1,1);
inf_core_core_q = zeros(length(quart)-1,1);
for j = 2:length(quart)
    inf_headline_q(j-1) = 400*log(cpi_headline_q(j)/cpi_headline_q(j-1));
    inf_core_q(j-1) = 400*log(cpi_core_q(j)/cpi_core_q(j-1));
    inf_boj_core_q(j-1) = 400*log(cpi_boj_core_q(j)/cpi_boj_core_q(j-1));
    inf_core_core_q(j-1) = 400*log(cpi_core_core_q(j)/cpi_core_core_q(j-1));
end

% Write to Excel
filename = 'AnnInf_q.xlsx';
inf_q = [quart(2:end)',inf_headline_q,inf_core_q,inf_boj_core_q,inf_core_core_q];
xlswrite(filename,inf_q)


% Plotting
inx = find(quart(2:end) == 1981);
inf_q_plot =inf_q(inx:end,:);

figure;
subplot(2,2,1)
box on
hold on
grid on
plot(inf_q_plot(:,1),inf_q_plot(:,3),'k','LineWidth',1.5)
xlim([inf_q_plot(1,1) inf_q_plot(end,1)])
ylim([-5 10])
title('Core')

subplot(2,2,2)
box on
hold on
grid on
plot(inf_q_plot(:,1),inf_q_plot(:,2),'k','LineWidth',1.5)
xlim([inf_q_plot(1,1) inf_q_plot(end,1)])
ylim([-5 10])
title('Headline')


subplot(2,2,3)
box on
hold on
grid on
plot(inf_q_plot(:,1),inf_q_plot(:,5),'k','LineWidth',1.5)
xlim([inf_q_plot(1,1) inf_q_plot(end,1)])
ylim([-5 10])
title('Core-Core')

subplot(2,2,4)
box on
hold on
grid on
plot(inf_q_plot(:,1),inf_q_plot(:,4),'k','LineWidth',1.5)
xlim([inf_q_plot(1,1) inf_q_plot(end,1)])
ylim([-5 10])
title('BoJ-Core')




    
    
