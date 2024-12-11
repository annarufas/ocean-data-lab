
% ======================================================================= %
%                                                                         %
%                       Chla climatologies from                           %
%               NASA Aqua-MODIS, NASA SeaWiFS and OC-CCI                  % 
%                                                                         %
% This script reads in chlorophyll concentration data from NASA's Ocean   % 
% Color website (downloaded manually following the link below) and from   %
% the Ocean Colour Climate Change Initiative (OC-CCI) (downloaded using   %
% the script "downloadChlaFromOCCI.m".                                    %
%                                                                         %
% NASA's Ocean Color website dataset characteristics:                     %
%   - https://oceancolor.gsfc.nasa.gov/l3/order/                          % 
%   - Time resolution: L3 monthly climatology                             %
%   - Space resolution: 4 km (8640 x 4320) (Aqua-MODIS),                  %
%                       9 km (4320 x 2160) (SeaWiFS)                      %
%   - Version: 2002–2024 (single sensor, Aqua-MODIS)                      %
%              1998-2010 (single sensor, SeaWiFS)                         %
%   - Units: mg m-3                                                       %
%                                                                         %
% Notice that if you click on the dataset link above, for the Aqua-MODIS  %
% sensor mapped monthly climatology, two months (July & August) are       % 
% missing from the collection of files. I accessed them following these   %
% instructions: https://forum.earthdata.nasa.gov/viewtopic.php?t=4784.    %
%                                                                         %
% OC-CCI dataset characteristics:                                         %
%   - https://www.oceancolour.org                                         % 
%   - Time resolution: monthly                                            %
%   - Space resolution: 4 km (8640 x 4320)                                %
%   - Version: v6, 1997–2024 (multi-sensor)                               %
%   - Units: mg m-3                                                       %
%                                                                         %
% This script uses this external function:                                %
%       processNASAsensorData - custom function                           %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 4 Nov 2024                                    %
%                                                                         %
% ======================================================================= %

% Clear workspace, close figures, and add paths to plotting resources
close all; clear all; clc
addpath(genpath(fullfile('code')));
addpath(genpath(fullfile('resources','external'))); 
addpath(genpath(fullfile('resources','internal'))); 
addpath(genpath(fullfile('figures')))

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

fullpathInputChlaOccciDir = fullfile('data','raw','chla','OCCCI');
fullpathInputChlaAquaModisDir = fullfile('data','raw','chla','AquaMODIS');
fullpathInputChlaSeawifsDir = fullfile('data','raw','chla','SeaWiFS');

fullpathOutputChlaAquaModisFile = fullfile('data','processed','chla_aquamodis.mat');
fullpathOutputChlaSeawifsFile = fullfile('data','processed','chla_seawifs.mat');
fullpathOutputChlaOccciFile = fullfile('data','processed','chla_occci.mat');

isDownloadOccciFilesCompleted = 1; % if 0, execute the function "downloadChlaFromOCCCI"

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - READ IN CHLA FILES FROM AQUA-MODIS
% -------------------------------------------------------------------------

% Get information
S = ncinfo(fullfile(fullpathInputChlaAquaModisDir,... 
    ['AQUA_MODIS.' '20020701_20230731.L3m.MC.CHL.chlor_a.4km' '.nc']));

% Read in the data
[chla,lat,lon] = processNASAsensorData(fullpathInputChlaAquaModisDir,...
    'AQUA_MODIS.','chlor_a','20020701_20230731.L3m.MC.CHL.chlor_a.4km');

% Check for spurious data points
figure(); histogram(chla(:), 100);

% Save the data
chla_lat = lat;
chla_lon = lon;
save(fullpathOutputChlaAquaModisFile,'chla','chla_lat','chla_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputChlaAquaModisFile,[],'mg m^{-3}',...
    0,1,true,'fig_monthly_chla_aquamodis','Chla Aqua-MODIS')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - READ IN CHLA FILES FROM SEAWIFS
% -------------------------------------------------------------------------

% Get information
S = ncinfo(fullfile(fullpathInputChlaSeawifsDir,... 
    ['SEASTAR_SEAWIFS_GAC.' '19970901_20100930.L3m.MC.CHL.chlor_a.9km' '.nc']));

% Read in the data
[chla,lat,lon] = processNASAsensorData(fullpathInputChlaSeawifsDir,...
    'SEASTAR_SEAWIFS_GAC.','chlor_a','19970901_20100930.L3m.MC.CHL.chlor_a.9km');

% Check for spurious data points
figure(); histogram(chla(:), 100);

% Save the data
chla_lat = lat;
chla_lon = lon;
save(fullpathOutputChlaSeawifsFile,'chla','chla_lat','chla_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputChlaSeawifsFile,[],'mg m^{-3}',...
    0,1,true,'fig_monthly_chla_seawifs','Chla SeaWiFS')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - READ IN CHLA FILES FROM OC-CCI
% -------------------------------------------------------------------------

if ~isDownloadOccciFilesCompleted
    downloadChlaFromOCCCI(fullpathInputChlaOccciDir)
end

% Get a list of all files in the directory
ncFiles = dir(fullfile(fullpathInputChlaOccciDir, '*.nc'));

% Get information
S = ncinfo(fullfile(ncFiles(1).folder, ncFiles(1).name));

% Read latite/longitude
lon = double(ncread(fullfile(ncFiles(1).folder, ncFiles(1).name),'lon')); 
lat = double(ncread(fullfile(ncFiles(1).folder, ncFiles(1).name),'lat'));

% Get file names grouped by month
groupedFiles = struct();
for iFile = 1:length(ncFiles)
    filename = ncFiles(iFile).name;
    filenameParts = split(filename, '-');
    yearStr = filenameParts{7}(1:4); % extract year (e.g., '1997')
    monthStr = filenameParts{7}(5:6); % extract month (e.g., '09')
    yearVal = str2double(yearStr);
    monthVal = str2num(monthStr); 
    if ~isfield(groupedFiles, ['month' num2str(monthVal)])
        groupedFiles.(['month' num2str(monthVal)]) = {};
    end
    groupedFiles.(['month' num2str(monthVal)]){end+1} = filename;
end

% For each month, read the data from the files and calculate the mean
chlorophylla = zeros(numel(lon),numel(lat),12); 

for iMonth = 1:12
    fileListForMonth = groupedFiles.(['month' num2str(iMonth)]);
    
    % Loop over files for the current month
    for i = 1:length(fileListForMonth)
        filename = fileListForMonth{i};
        chlorData = ncread(fullfile(fullpathInputChlaOccciDir,filename), 'chlor_a'); 
        
        % On the first iteration, initialise the sum and count
        if i == 1
            sumData = zeros(size(chlorData)); 
            countData = zeros(size(chlorData));
        end
        
        % Operate only for non-NaN values
        validMask = ~isnan(chlorData);                                  % logical mask of valid values
        sumData(validMask) = sumData(validMask) + chlorData(validMask); % sum only valid data
        countData(validMask) = countData(validMask) + 1;                % increment count only for valid data
        
    end 
    chlorophylla(:,:,iMonth) = sumData ./ countData; % mean data   
end

% Some arrangements
chlorophylla(chlorophylla < 0) = 0; % we cannot have negative values

% Sort latitudes to have monotonically increasing values
[lat_sort, sortIdx] = sort(lat);
chla_sort = chlorophylla(:,sortIdx,:);

% Swap lon and lat dimensions to get lat x lon x depth
chla_sort_perm = permute(chla_sort, [2, 1, 3]);

% Check for spurious data points
figure(); histogram(chla_sort_perm(:),100);

% Save the data
chla = chla_sort_perm;
chla_lat = lat_sort;
chla_lon = lon;
save(fullpathOutputChlaOccciFile,'chla','chla_lat','chla_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputChlaOccciFile,[],'mg m^{-3}',...
    0,1,true,'fig_monthly_chla_occci','Chla OC-CCI')
