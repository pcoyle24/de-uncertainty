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
starttime = 1998.0;
endtime = 2023.5;
time = (starttime:0.25:endtime)';

%% Inflation Data
sheet_names = {'am01-1'};
values = {'I15:N662'}; % from 1970M1 to 2023M12
file_name = strcat(datapath,'raw/cpi.xlsx');
cpi_data = xlsread(file_name,char(sheet_names(1)),char(values));

% Rows are:
% 1. All Items
% 2. All Items (less fresh food)
% 3. All Items (less fresh food and energy)
% 4. All Items (Less fresh food -- less alcoholic bev -- and energy)

% Not Seasonally Adjusted
% cpi_data = [cpi_data(:,1),...
%             cpi_data(:,2),...
%             cpi_data(:,5), ...
%             cpi_data(:,6)];

% Seasonally Adjust using X13
cpi_data = [X13SA(cpi_data(:,1), 'M', '1970M1'),...
            X13SA(cpi_data(:,2), 'M', '1970M1'),...
            X13SA(cpi_data(:,5), 'M', '1970M1'), ...
            X13SA(cpi_data(:,6), 'M', '1970M1')];

time_cpi = (1970:0.25:2023.75)';

cpi_headline = MonthToQuarter(cpi_data(:,1), time_cpi(1), time_cpi(end), 'Avg');
cpi_core = MonthToQuarter(cpi_data(:,2), time_cpi(1), time_cpi(end), 'Avg');
cpi_boj_core = MonthToQuarter(cpi_data(:,3), time_cpi(1), time_cpi(end), 'Avg');
cpi_core_core = MonthToQuarter(cpi_data(:,4), time_cpi(1), time_cpi(end), 'Avg');

time_inf = (1971:0.25:2023.75)';
inf_sinx = find(time_inf == starttime);
inf_einx = find(time_inf == endtime);
% inf = [100*log(cpi_headline(5:end)./cpi_headline(1:end-4)), ...
%     100*log(cpi_core(5:end)./cpi_core(1:end-4)), ...
%     100*log(cpi_boj_core(5:end)./cpi_boj_core(1:end-4)), ...
%     100*log(cpi_core_core(5:end)./cpi_core_core(1:end-4))];

inf = [400*log(cpi_headline(2:end)./cpi_headline(1:end-1)), ...
    400*log(cpi_core(2:end)./cpi_core(1:end-1)), ...
    400*log(cpi_boj_core(2:end)./cpi_boj_core(1:end-1)), ...
    400*log(cpi_core_core(2:end)./cpi_core_core(1:end-1))];

%% Output Gap Data (BOJ)
sheet_names = {'Chart'};
values = {'B6:B168'}; % from 1983Q1 to 2023Q3
file_name = strcat(datapath,'raw/gap_boj.xlsx');
gap_boj = xlsread(file_name,char(sheet_names(1)),char(values));
time_gap_boj = (1983:0.25:2023.50)';
gap_boj_sinx = find(time_gap_boj == starttime);
gap_boj_einx = find(time_gap_boj == endtime);

%% Output Gap Data (Cabinet)
sheet_names = {'Chart'};
values = {'C7:C181'}; % from 1980Q1 to 2023Q3
file_name = strcat(datapath,'raw/gap_cabinet.xlsx');
gap_cabinet = xlsread(file_name,char(sheet_names(1)),char(values));
time_gap_cabinet = (1980:0.25:2023.50)';
gap_cabinet_sinx = find(time_gap_cabinet == starttime);
gap_cabinet_einx = find(time_gap_cabinet == endtime);

%% GDP Deflator
values = {'B8:B127'}; % from 1994Q1 to 2023Q4
file_name = strcat(datapath,'raw/defl.xlsx');
defl_data = xlsread(file_name,char(values)); %Already Seasonally Adjusted

% defl = 100*log(defl_data(5:end)./defl_data(1:end-4));
defl = 400*log(defl_data(2:end)./defl_data(1:end-1));

time_defl = (1995:0.25:2023.75)';
defl_sinx = find(time_defl == starttime);
defl_einx = find(time_defl == endtime);


%% Nominal Interest Rate
values = {'C70:C531'}; % from 1985M7 to 2023M12
file_name = strcat(datapath,'raw/intrate.xlsx');
int_data = xlsread(file_name,char(values));
time_int = (1985.5:0.25:2023.75)';
int_sinx = find(time_int == starttime);
int_einx = find(time_int == endtime);

int = MonthToQuarter(int_data, time_int(1), time_int(end), 'Avg');


%% Final Data
data = [time, ...
        inf(inf_sinx:inf_einx,:),...
        gap_boj(gap_boj_sinx:gap_boj_einx,:)/4,...
        gap_cabinet(gap_cabinet_sinx:gap_cabinet_einx,:)/4,...
        defl(defl_sinx:defl_einx,:), ...
        int(int_sinx:int_einx,:)];

T = table(data(:,1),...
      data(:,2),...
      data(:,3),...
      data(:,4),...
      data(:,5),...
      data(:,6),...
      data(:,7),...
      data(:,8),...
      data(:,9),...
      'VariableNames', ["Time",...
                        "Headline Inflation",...
                        "Core Inflation",...
                        "Core Inflation (BOJ)",...
                        "Core Core Inflation",...
                        "Gap (BOJ)",...
                        "Gap (Cabinet)",...
                        "GDP Defl",...
                        "Int Rate"]);

file_name = strcat(datapath,'cleaned/data.xlsx');
command = strcat('rm -r', " ", file_name);
system(command);
writetable(T,file_name,'Sheet',1,'Range','A1')


