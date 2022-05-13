function varargout = midpoint(x)
    arguments (Repeating)
        x (1,:)
    end
% find the midpoint of indicies in varargin
    varargout = cell(length(x),1);
    for ii = 1:length(x)
        varargout{ii} = mean(x{ii});
    end


end