% dtick(x) - specifies datelabels in format 'x' for x-axis in 5 ticks
% fbl 03/08/05

function dtick(varargin)

if ~isempty(varargin)
    x = cell2mat(varargin(1));
else
    x = 24;      % date format %24    'mmm.dd.yyyy HH:MM'
end  

lim = get(gca,'xlim');
datestr(xlim);
set(gca,'xtick',[lim(1):(lim(2)-lim(1))/4:lim(2)]);
datetick('x',x,'keepticks');

end


% 0     'dd-mmm-yyyy HH:MM:SS'
% 15    'HH:MM'
% 19    'dd/mm'
% 21    'mmm.dd.yyyy HH:MM:SS'
