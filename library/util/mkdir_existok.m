function mkdir_existok(path)

if ~isfolder(path)
    mkdir(path)
end