function sa_data = X13SA(Data, Freq, StartPeriod)
    % Tutorial: See https://vermandel.fr/2022/12/30/deseasonalize-time-series-with-x13-in-matlab/
    % Must have Dynare installed which provides the dseries function and
    % x13 function. 

    % obtain Dynare directory
    get_dynare_src = strrep(which('dynare'),'dynare.m','');
    
    % load dseries object
    if ispc
        addpath([get_dynare_src 'modules\dseries\src\'],[get_dynare_src 'missing\rows_columns\'])
    else
        addpath([get_dynare_src '/modules/dseries/src/'],[get_dynare_src '/missing/rows_columns/'])
    end
    initialize_dseries_class();
    
    % convert data into dseries object
    ts = dseries(Data,StartPeriod);
    
    % create the x13 object
    o = x13(ts);
    % adjust options
    % o.transform('function','log');
    if strcmp(Freq, 'Q')
        o.arima('model',' (0 1 1)4');
    elseif strcmp(Freq, 'M')
        o.arima('model',' (0 1 1)12');
    end
        
    o.x11('save','(d10)');
    % run
    o.run();

    % delete useless files program produces
    hk = char('d10','err','log','out','spc');
    for i0 = 1:size(hk,1)
        eval(['delete ' o.results.name '.' hk(i0,:) ';']);
    end
    
    % extract the multiplicative seasonal pattern
    season_y = o.results.d10;
    
    sa_data = (o.y.data)./(season_y.data);
end

