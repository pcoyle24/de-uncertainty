function []  = axeslim(handle,ybuffer,digits,varargin)

% axeslim adjusts plot axes for figure or subplot
% []  = axeslim(handle,ybuffer,digits,varargin)
% Inputs:
%   handle     :   Figure handle
%   ybuffer     :   Determines the buffer on the yaxis
%   digits      :   Rounds ylimits to the specified decimal
%   varargin       
%    xparts     :   number of xticks  
%    precision  :   rounds yticks to the specified decile

% Authors:
%   Alexander Richter (richtera@indiana.edu)
%   Nathaniel Throckmorton (nathrock@indiana.edu)
%   Indiana University - Department of Economics

% Get axis handles
if strcmp(get(handle,'Type'),'figure')
    children = get(handle,'Children');
    idx = 1;
    for i = 1:length(children)
        if strcmp(get(children(i),'Type'),'axes')
            axs(idx) = children(i);
            idx = idx + 1;
        end
    end
elseif strcmp(get(handle,'Type'),'axes')
    axs = handle;
end
% Set axis properties
for i = 1:length(axs)
    % Get line handles
    children = get(axs(i),'Children');
    idx = 1;
    for j = 1:length(children)
        if strcmp(get(children(j),'Type'),'line')
            lines(idx) = children(j);
            idx = idx + 1;
        end
    end
    % Grab data for each line
    ydata = [];
    xdata = [];    
    for j = 1:length(lines)        
        ydata = [ydata get(lines(j),'Ydata')];
        xdata = [xdata get(lines(j),'Xdata')];
    end
    % Determine axis limits
    ylim1 = round(10^digits*min(ydata))/10^digits;  
    ylim2 = round(10^digits*max(ydata))/10^digits;
    step = round(10^digits*abs(ylim1-ylim2)*ybuffer)/10^digits;
    if ylim1 < 0
        ylim1 = ylim1 - step;
    else
        ylim1 = max(0,ylim1 - step);
    end
    ylim2 = ylim2 + step;
    xlim1 = ceil(100*min(xdata))/100;
    xlim2 = floor(100*max(xdata))/100;
    % Set axis limits
    if abs(ylim1-ylim2) > 1e-4
        axis(axs(i),[xlim1,xlim2,ylim1,ylim2]);
    end
    if ~isempty(varargin)
        xparts = varargin{1};            
        % Set xaxis ticks
        xstep = round(100*(xlim2-xlim1)./xparts)/100;
        set(axs(i),'xtick',xlim1:xstep:xlim2);
        set(axs(i),'xticklabel',xlim1:xstep:xlim2);  
        set(axs(i),'xlim',[xlim1,xlim2])
        % Remove scientific notation from axis
        if numel(varargin)==2
            precision = varargin{2};
            tick = get(axs(i),'ytick');       
            set(axs(i),'yticklabel',round(tick*precision)/precision,'yticklabelmode','manual')         
            set(gca,'YTickLabelMode','auto')
            set(gca,'YTickMode','auto')
        end
    end
end