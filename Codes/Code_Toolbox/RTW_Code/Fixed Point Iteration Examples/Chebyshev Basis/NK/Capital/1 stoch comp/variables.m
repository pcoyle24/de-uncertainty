function V = variables

% V = variables
%   Names variables and assigns locations
% Output:
%   V : Structure of variable names and locations

% Variables names
V.names = { 'n'          %1  Labor hours
            'w'          %2  Real wage rate
            'c'          %3  Consumption
            'm'          %4  Real money balances
            'r'          %5  Net nominal interest rate 
            'pi'         %6  Gross inflation
            'y'          %7  Output
            'psi'        %8  Marginal cost
            'tau'        %9  Tax rate
            'b'          %10 Real Debt
            'rk'         %11 Rental rate of capital
            'k'          %12 Capital
            'i'          %13 Investment
            'a'          %14 Government Liabilities 
            'tr'         %15 Tax Revenue
          };

% Variables titles       
V.desc = {  'Labor Hours (%)'
            'Real Wage (%)'
            'Consumption (%)'
            'Real Money Balances (%)'
            'Nominal Interest Rate (% point)'
            'Inflation (% point)'
            'Output (%)'
            'Marginal Cost (%)'
            'Lump-Sum Tax (% point)'
            'Real Debt (%)'
            'Rental Rate of Capital (%)'
            'Capital (%)'
            'Investment (%)'
            'Government Liabilities (%)'
            'Tax Revenues (%)'
         };

%Shocktypes: 
V.shocktypes = {'fp'};  
% Forecast errors
V.foretypes = {'c','pi','rk','tau'}; 
             
% Number of variables
V.nvar = length(V.names);
% Number of shocks
V.nshock = length(V.shocktypes);
% Number of forecast errors
V.nfore = length(V.foretypes);

% Establish variable index
for j = 1:size(V.names,1)
   eval(['V.' V.names{j} ' = j;']);      
end
% Establish shock index
for j = 1:V.nshock
    eval(['V.eps' V.shocktypes{j} ' = j;'])
end
% Establish forecast error index
for j = 1:V.nfore
    eval(['V.fe' V.foretypes{j} ' = j;'])
end