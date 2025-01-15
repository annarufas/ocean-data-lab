
% ======================================================================= %
%                                                                         %
%                  SST climatology from AVHRR Pathfinder                  % 
%                                                                         %
% This script reads in sea surface temperature (SST) data from AVHRR      %
% Pathfinder v5.0.                                                        %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:AVHRR_Pathfinder-NODC-v5.0-climatologies
%   - Time resolution: monthly climatology                                %
%   - Space resolution: 4 km (4096 x 8192)                                %
%   - Version: v5, 1985–2001                                              %
%   - Units: ºC                                                           %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
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

filenameOutputSst = 'sst_pathfinder_v5.mat';
fullpathInputSstDir = fullfile('data','raw','SST','AVHRR_Pathfinder_v5');
fullpathOutputSstFile = fullfile('data','processed',filenameOutputSst);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

% The data comes in unit16 format, to convert to double apply the
% scale factor (slope) of 0.075 and offset (y-intercept) of -3, as instructed
% in the attributes of the metadata and here: https://grasswiki.osgeo.org/wiki/AVHRR

SCALE_FACTOR = 0.075;
OFFSET_FACTOR = -3;
LAND_MASK = OFFSET_FACTOR + SCALE_FACTOR;

S = hdfinfo(fullfile(fullpathInputSstDir,'month01_combined.hdf'));

% Read in data
lon = hdfread(fullfile(fullpathInputSstDir,...
    ['month' sprintf('%02d',1) '_combined.hdf']),'Longitude')'; 
lat = hdfread(fullfile(fullpathInputSstDir,...
    ['month' sprintf('%02d',1) '_combined.hdf']),'Latitude')'; 

seasurfacetemp = zeros(numel(lat),numel(lon),12);
for iMonth = 1:12 
    values = double(hdfread(fullfile(fullpathInputSstDir,...
        ['month' sprintf('%02d',iMonth) '_combined.hdf']),'Clim_SST_Filled'));
    seasurfacetemp(:,:,iMonth) = SCALE_FACTOR.*values + OFFSET_FACTOR;
end

% Some arrangements (mask land)
seasurfacetemp(seasurfacetemp <= LAND_MASK) = NaN;

% Sort latitude to have monotonically increasing values
[lat_sort, sortIdx] = sort(lat);
sst_sort = seasurfacetemp(sortIdx,:,:);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CHECK VALUES AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(sst_sort(:),100); ylim([0 3e7])

% Save the data
sst = sst_sort;
sst_lon = lon;
sst_lat = lat_sort;
save(fullpathOutputSstFile,'sst','sst_lat','sst_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputSstFile,[],'ºC',...
    -2,35,true,'fig_monthly_sst_pathfinder','SST Pathfinder')
