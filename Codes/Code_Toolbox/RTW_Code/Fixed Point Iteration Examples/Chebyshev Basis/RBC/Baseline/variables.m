function V = variables

% V = variables
%   Names variables and assigns locations
% Output:
%   V : Structure of variable names and locations

% Variables names
V.names = { 'n'          %1  Labor hours
            'w'          %2  Real wage rate
            'c'          %3  Consumption
            'y'          %4  Output
            'rk'         %5  Rental rate of capital
            'k'          %6  Capital
            'i'          %7  Investment
            'z'          %8  Total Factor Productivity
          };

% Variables titles        
V.desc = {  'Labor Hours'
            'Real Wage Rate'
            'Consumption'
            'Output'
            'Rental Rate of Capital'
            'Capital'
            'Investment'
            'Productivity'
         };

%Shocktypes:
V.shocktypes = {'z'};  
% Forecast errors
V.foretypes = {'c','rk'}; 
             
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