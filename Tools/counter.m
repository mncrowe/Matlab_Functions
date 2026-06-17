function counter(t, p, units, txt)
% Gives current progress and time remaining
%
% - t: current time, e.g. determined using toc function
% - p: current progress (normalised to 1)
% - units: 0 - determine automatically, 1 - secs, 2 - mins, 3 - hours, 4 - days (default: 0)
% - txt: optional additional text to display after remaining time

arguments
    t (1,1) double
    p (1,1) double
    units (1,1) int64 = 0
    txt char = '';
end

if units > 4; units = 0; end
txt = [', ' txt];

t_rem = t / p * (1-p);

if units == 0
    units = 1;
    if t_rem > 60; units = 2; end
    if t_rem > 60 * 60; units = 3; end
    if t_rem > 60 * 60 * 24; units = 4; end
end

switch units
    case 1
        units_str = 'secs';
    case 2
        units_str = 'mins';
	    t_rem = t_rem / 60;
    case 3
        units_str = 'hours';
	    t_rem = t_rem / (60 * 60);
    otherwise
        units_str = 'days';
	    t_rem = t_rem / (60 * 60 * 24);
end

t_rem=round(100 * t_rem) / 100;

disp([num2str(round(p * 1000) / 10) ' % complete, remaining time ~ ' num2str(t_rem) ' ' units_str txt]);

end
