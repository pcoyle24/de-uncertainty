function V = variables

% V = variables
%   Names variables and assigns locations
% Output:
%   V : Structure of variable names and locations

% Variables names
V.names = { 'c'     %(1) Consumption
            'r'     %(2) Gross Nominal Interest Rate
            'inf'   %(3) Gross Inflation
            'n'     %(4) Labor
            'del'   %(5) Demand Shock
          };

% Variables titles       
V.desc = {  'Consumption (%)'
            'Gross Nominal Interest Rate (%)'
            'Gross Inflation (%)'
            'Labor (%)'
            'Demand Shock (% point)'
         };

%Shocktypes: (demand)
V.shocktypes = {'d'};  
% Forecast errors %What we are solving for in the model
V.foretypes = {'c','n','inf'}; 
             
% Number of variables
V.nvar = length(V.names);
% Number of shocks
V.nshock = length(V.shocktypes);
% Number of forecast errors
V.nfore = length(V.foretypes);

% Establish variable index
for j = 1:size(V.names,1)
%for j = 1:length(V.names) %Try this - it might work?
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