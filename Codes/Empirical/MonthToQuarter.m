function out = MonthToQuarter(data, startq, endq, method)
    
    quart = startq:.25:endq;
    out = zeros(length(quart),1);

    % Convert Monthly Series to Quarterly Series
    j = 1;
    for i = 1:3:length(data)-1
        if strcmp(method, 'Avg')
            out(j) = sum(data(i:i+2))/3;
        elseif strcmp(method, 'EndOfQtr')
            out(j) = data(i+2);
        end
        j = j+1;
    end

end
