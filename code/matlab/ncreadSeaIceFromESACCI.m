
% ======================================================================= %
%                                                                         %
%           Sea ice concentration climatology from ESA-CCI                %
%                                                                         %
% This script reads in sea ice concentration data from the European Space %
% Agency Climate Change Initiative (ESA-CCI).                             %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://catalogue.ceda.ac.uk/uuid/5f75fcb0c58740d99b07953797bc041e  % 
%   - Time resolution: month                                              %
%   - Space resolution: 25 km resolution (N. and S. Hemisphere)           %
%   - Version: v2.1, 2002-2017                                            %
%   - Units: percentage                                                   %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 11 Jan 2025                                   %
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

filenameOutputSeaIce = 'icefrac_esacci.mat';
fullpathInputSeaIceDir = fullfile('data','raw','seaice','ESACCI');
fullpathOutputSeaIceFile = fullfile('data','processed',filenameOutputSeaIce);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

% Get a list of all files in the directory
ncFiles = dir(fullfile(fullpathInputSeaIceDir)); 
ncFilesNH = ncFiles(contains({ncFiles.name}, 'NH'));
ncFilesSH = ncFiles(contains({ncFiles.name}, 'SH'));

% Initialise output structure
Dout = struct();

for i = 1:2
    % Select file group
    if i == 1
        fileList = ncFilesNH;
        region = 'NH'; % Northern Hemisphere
    elseif i == 2
        fileList = ncFilesSH;
        region = 'SH'; % Southern Hemisphere
    end

    % Extract months from filenames
    fileMonths = zeros(length(fileList), 1);
    for iFile = 1:length(fileList)
        % Split filename and extract the date part (assumes NH-YYYYMMDD format)
        fileNameParts = split(fileList(iFile).name, '-');
        if length(fileNameParts) > 6
            datePart = fileNameParts{7}; % e.g., "20020603"
            monthStr = datePart(5:6); % extract "MM" part
            fileMonths(iFile) = str2double(monthStr); 
        end
    end

    % Sort files by month
    [~, sortedIdx] = sort(fileMonths);
    fileList = fileList(sortedIdx); % reorder file list
    
    % Read lon/lat data from an example file
    exampleFile = fullfile(fileList(1).folder, fileList(1).name);
%     S = ncinfo(exampleFile); % short summary
    lon = double(ncread(exampleFile, 'lon'));
    lat = double(ncread(exampleFile, 'lat'));

    % Read data from all files
    D = zeros([size(lon) length(fileList)]); % lat and lon have the same size
    for iFile = 1:length(fileList)
        filePath = fullfile(fileList(iFile).folder, fileList(iFile).name);
        D(:,:,iFile) = double(ncread(filePath, 'ice_conc')); % units are "%"
    end
    
    % Initialise array to store monthly averages
    Dmonthly = zeros([size(lon), 12]);
    monthCounts = zeros(12, 1); % to count number of files per month
    
    % Compute monthly averages
    for iFile = 1:length(fileList)
        currentMonth = fileMonths(iFile);
        Dmonthly(:,:,currentMonth) = Dmonthly(:,:,currentMonth) + D(:,:,iFile);
        monthCounts(currentMonth) = monthCounts(currentMonth) + 1;
    end
    
    % Divide by counts to get averages
    for iMonth = 1:12
        Dmonthly(:,:,iMonth) = Dmonthly(:,:,iMonth) / monthCounts(iMonth);
    end
    
    % Store data in the output structure
    Dout.(region).data = Dmonthly;
    Dout.(region).lat = lat;
    Dout.(region).lon = lon;
    
end
        
% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CONCATENATE TO CREATE GLOBAL ARRAY
% -------------------------------------------------------------------------

% Concatenate the Southern Ocean and Northern Hemisphere data to form a global map
globalArray = [Dout.('SH').data; Dout.('NH').data]; 
globalLat   = [Dout.('SH').lat;  Dout.('NH').lat]; 
globalLon   = [Dout.('SH').lon;  Dout.('NH').lon];

globalArray(globalArray == 0) = NaN;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CHECK VALUES AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
% figure(); histogram(globalArray(:),100);

% Save the data
icefrac_lon = globalLon;
icefrac_lat = globalLat;
icefrac = globalArray;
save(fullpathOutputSeaIceFile,'icefrac','icefrac_lat','icefrac_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputSeaIceFile,[],'%',...
    0,100,true,'fig_monthly_icefrac_esacci','Ice fraction ESA-CCI')
