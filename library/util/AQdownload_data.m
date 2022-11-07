function [dataset] = AQdownload_data(dataname)
% Get an AirQuant example dataset
%
% Download an AirQuant example dataset if it hasn't already been downloaded
%
% .. todo:
%   * Automatically infer download size from HTTP head request of
%       link.
%   * Add exception to download whole file if interrupted.
%
% .. note:
%   Example datasets are saved in the default AirQuant directory
%     even if specified otherwise by the user.
%
%
% Args:
%   dataname('char'): name of dataset to download.
%
% Example:
%   >>> AQdownload_data('chestct')
%

% check if its already downloaded
AQroot = AirQuantAddPath();
dataset = 'airquant';
datadir = fullfile(AQroot,'data', dataset);
if ~exist(datadir,'dir')
    mkdir(datadir)
end
filename_noext = fullfile(datadir,dataname);
filename_withext = [filename_noext, '.tar.gz'];
if ~isempty(dir([filename_noext,'*']))
    disp(['AirQuant Dataset: "',dataname, '" already exists.'])
    return
end

% get url and download
switch dataname
    case 'chestct'
        url = 'https://dl.dropboxusercontent.com/s/g94zhgyvknn512w/chestct.tar.gz?dl=0';
        size = 257;
end

disp(['Downloading AirQuant dataset, "', dataname, '" of ', num2str(size), ' MB. This may ',...
    'take a while depending on your connection.'])
websave(filename_withext, url);
disp('This case is from an opensource dataset, please see AirQuant readme for credits.')

% unzip and delete compressed
untar(filename_withext,datadir)
delete(filename_withext)

end
