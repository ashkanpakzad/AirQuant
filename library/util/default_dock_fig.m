function default_dock_fig(flag)
% set default behaviour for new figures to be docked.
if nargin < 1
    flag = 1;
end

if flag == 1
    set(0,'DefaultFigureWindowStyle','docked')
else
    set(0,'DefaultFigureWindowStyle','normal')
end

end