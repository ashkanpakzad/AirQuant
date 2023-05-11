function fig_save(fig_handle,path,ext)
% saves figure as both .mat and .png (determined by ext).
%
% path must not have an ext.
%

if nargin < 3
    ext = '.png';
end

% save as MATLAB proprietary .fig
savefig(fig_handle, path)
% save as image file
exportgraphics(fig_handle, strcat(path,ext));
end

