function [D_sort_perm,lat_sort,lon] = processNASAsensorData(inputDir,prefixTag,...
    varName,exampleFile)

% PROCESSNASASENSORDATA This function function that encapsulates the 
% flow to extract data from NASA's Ocean Color sensors.
% 
%   INPUT: 
%       inputDir    - full path to data input directory
%       prefixTag   - prefix common to all climatology files
%       varName     - name of the avriable we want to extract
%       exampleFile - example file tor read in latitudes and longitudes
%                     (common arrangament to all files in the climatology)
%
%   OUTPUT:
%       D_sort_perm - variable
%       lat_sort    - latitudes vector arranged -90 to 90
%       lon         - longitudes vector arranged -180 to 180
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 4 Nov 2024 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

% Get a list of all files in the directory
ncFiles = dir(fullfile(inputDir, [prefixTag '*.nc']));

% Sort files based on month extracted from the filename
fileMonths = zeros(length(ncFiles),1);
for iMonth = 1:length(ncFiles)
    fileNameParts = split(ncFiles(iMonth).name, '.');
    startMonthStr = fileNameParts{2}(5:6); % extract start month (e.g., '06')
    fileMonths(iMonth) = str2num(startMonthStr); % convert date string to serial date number
end
[~, sortedIdx] = sort(fileMonths); % sort the files based on the extracted dates
ncFiles = ncFiles(sortedIdx); % reorder ncFiles based on sorted indices

% Read in the data
lon = double(ncread(fullfile(inputDir,strcat(prefixTag,exampleFile,'.nc')), 'lon'));
lat = double(ncread(fullfile(inputDir,strcat(prefixTag,exampleFile,'.nc')), 'lat'));

D = zeros(numel(lon),numel(lat),length(ncFiles)); 
for iMonth = 1:length(ncFiles)
    filePath = fullfile(ncFiles(iMonth).folder,ncFiles(iMonth).name);
    D(:,:,iMonth) = double(ncread(filePath,varName));
end

% Set negative values to zero
D(D < 0) = 0; 

% Sort latitudes to have monotonically increasing values
[lat_sort,sortIdx] = sort(lat);
D_sort = D(:,sortIdx,:);

% Swap lon and lat dimensions to get lat x lon x time
D_sort_perm = permute(D_sort, [2, 1, 3]);

end % processNASAsensorData

