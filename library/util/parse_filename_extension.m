function filename = parse_filename_extension(filename, requiredext)
    % check extension is requiredext, if not change/set extenesion.
    [filepath,name,ext] = fileparts(filename);
    if ~strcmp(ext, requiredext)
        filename =  fullfile(filepath, strcat(name,requiredext));
    end
end