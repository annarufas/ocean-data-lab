
% ======================================================================= %
%                                                                         %
%                     PAR0 climatology from GlobColour                    % 
%                                                                         %
% This script reads in photosynthetic available radiation at the surface  %
% ocean (PAR0) from the European Space Agency (ESA) GlobColour project    %
% (https://hermes.acri.fr), which delivers merged satellite-derived       %
% data products on ocean colour and related variables. Some merged L3     % 
% products are available through the Copernicus Marine Service (CMEMS)    % 
% Data Store, but not PAR0. Unlike individual sensor products on the      %
% Ocean Biology Processing Group (OBPG) website                           %
% (https://oceancolor.gsfc.nasa.gov), gaps due to the orbit and sensing   %
% geometry are reduced. PAR0 data was downloaded using the script         %
% "downloadOCfromGlobColour.ipynb".                                       %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://hermes.acri.fr/index.php?class=archive                      % 
%   - Time resolution: L3 monthly                                         %
%   - Space resolution: 100 km (360 x 180)                                %  
%   - Version 2024                                                        %
%         1997-2002 (single sensor, SeaWiFS)                              %
%         2002-2010 (merged sensor, SeaWiFS and Aqua-MODIS)               %
%         2010-2012 (single sensor, Aqua-MODIS)                           %
%         2012-2024 (merged sensor, Aqua-MODIS and VIIRS)                 %
%   - Units: einstein m-2 d-1 = mol phot. m-2 d-1 (input), W m-2 (output) %
%                                                                         % 
% This script uses this external function:                                %
%       downloadOCfromGlobColour.ipynb - custom function                  %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 22 Dec 2024                                   %
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

fullpathInputPar0Dir = fullfile('data','raw','GlobColour');
fullpathOutputPar0File = fullfile('data','processed','par0_globcolour.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

% Get a list of all files in the directory
ncFiles = dir(fullfile(fullpathInputPar0Dir, '*.nc'));

% Get information
S = ncinfo(fullfile(ncFiles(1).folder, ncFiles(1).name));

% Read latite/longitude
lon = double(ncread(fullfile(ncFiles(1).folder, ncFiles(1).name),'lon')); 
lat = double(ncread(fullfile(ncFiles(1).folder, ncFiles(1).name),'lat'));

% Get file names grouped by month
groupedFiles = struct();
for iFile = 1:length(ncFiles)
    filename = ncFiles(iFile).name;
    filenameParts = split(filename, '_');
    %yearStr = filenameParts{2}(1:4); % extract year (e.g., '1997')
    monthStr = filenameParts{2}(5:6); % extract month (e.g., '09')
    monthVal = str2num(monthStr); 
    if ~isfield(groupedFiles, ['month' num2str(monthVal)])
        groupedFiles.(['month' num2str(monthVal)]) = {};
    end
    groupedFiles.(['month' num2str(monthVal)]){end+1} = filename;
end

% For each month, read the data from the files and calculate the mean
photosyntheticactiveradiation = NaN(numel(lon),numel(lat),12); 

for iMonth = 1:12
    fileListForMonth = groupedFiles.(['month' num2str(iMonth)]);
    
    % Loop over files for the current month
    for i = 1:length(fileListForMonth)
        filename = fileListForMonth{i};
        data = ncread(fullfile(fullpathInputPar0Dir,filename), 'PAR_mean'); 
        
        % On the first iteration, initialise the sum and count
        if i == 1
            sumData = zeros(size(data)); 
            countData = zeros(size(data));
        end
        
        % Operate only for non-NaN values
        validMask = ~isnan(data) & data ~= -999;                   % logical mask of valid values
        sumData(validMask) = sumData(validMask) + data(validMask); % sum only valid data
        countData(validMask) = countData(validMask) + 1;           % increment count only for valid data
        
    end 
    photosyntheticactiveradiation(:,:,iMonth) = sumData ./ countData; % mean data
end

% Some arrangements
photosyntheticactiveradiation(photosyntheticactiveradiation < 0) = 0; % we cannot have negative values

% Sort latitudes to have monotonically increasing values
[lat_sort, sortIdx] = sort(lat);
par0_sort = photosyntheticactiveradiation(:,sortIdx,:);

% Swap lon and lat dimensions to get lat x lon x time
par0_sort_perm = permute(par0_sort, [2, 1, 3]);

% Units conversion: einstein m-2 d-1 = mol photons m-2 d-1 -> umol photons m-2 s-1 -> W m-2
par0_sort_perm_Wm2 = par0_sort_perm .* (1e6 ./ (3600 .* 24)) .* (3.90e-19 .* 6.02e23 ./ 1e6); 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - INSPECT AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(par0_sort_perm_Wm2(:),100);

% Save the data
par0 = par0_sort_perm_Wm2;
par0_lat = lat_sort;
par0_lon = lon;
save(fullpathOutputPar0File,'par0','par0_lat','par0_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputPar0File,[],'W m^{-2}',...
    0,200,true,'fig_monthly_par0_globcolour','PAR0 GlobColour')
