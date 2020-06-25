function it = itinfo(istart,tstart,dec,it,dist_max)

% [Tlast UTh UTm UTs] = itinfo(istart,tstart)
% 	Calculates algorithm and iteration duration
% Output:
%   it     : Iteration counter
% Inputs:
%   istart : Start of iteration (from tic)
%   tstart : Start of algorithm (from tic)
%   dec    : Number of decimal places on ss (seconds)

dec = 10^dec;
T = toc(tstart);
itss = round(toc(istart)*dec)/dec;
hh = floor(T/3600);
mm = floor(T/60)-60*hh;
ss = round((T-60*mm-3600*hh)*dec)/dec;

% Convert to strings
dist_max = num2str(dist_max);
itss = num2str(itss);
hh = num2str(hh);
mm = num2str(mm);
ss = num2str(ss);

display(['iter: ' num2str(it) ' dist: ' dist_max ...
        ' Tlast: ' itss 's TTotal: ' hh 'h' mm 'm' ss 's']);

% Increase iteration counter        
it = it + 1;          