
% ======================================================================= %
%                                                                         %
%           Cloud cover climatology from Pincus et al. (2008)             %
%                                                                         %
% This script reads in cloud cover data from the publication of Pincus et %
% al. (2008). The cliamtology was created using CMIP3 data. As of October %
% 2024, the link to the dataset referenced in the publication is broken.  % 
%                                                                         %
% Dataset characteristics:                                                %
%   - https://doi.org/10.1029/2007JD009334 (downloaded Aug 2019)          %                          
%   - Time resolution: monthly climatology                                %
%   - Space resolution: 2.5 × 2.5º (144 x 72)                             %
%   - Version: 2008                                                       %
%   - Units: 0-1 (input), oktas (output)                                  %            
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
%                                                                         %
% ======================================================================= %

% Clear workspace, close figures, and add paths to plotting resources
addpath(genpath(fullfile('code')));
addpath(genpath(fullfile('resources','external'))); 
addpath(genpath(fullfile('resources','internal'))); 
addpath(genpath(fullfile('figures')))

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

filenameInputCloud = 'obs_data_set_secondary.nc';
fullpathInputCloudFile = fullfile('data','raw','CloudCover','Pincus',filenameInputCloud);
filenameOutputCloudClimatology = 'cloudcover_pincus.mat';
fullpathOutputCloudFile = fullfile('data','processed',filenameOutputCloudClimatology);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

S = ncinfo(fullfile(fullpathInputCloudFile));

% Read in the data
lon = double(ncread(fullfile(fullpathInputCloudFile),'lon')); 
lat = double(ncread(fullfile(fullpathInputCloudFile),'lat'));
cloudfrac = double(ncread(fullfile(fullpathInputCloudFile),'clt'));

% Some arrangements
cloudfrac(cloudfrac < 0) = 0; % we cannot have negative values

% Unit conversion, where 0 oktas = no clouds, 8 oktas = full cloud cover
cloudoktas = cloudfrac.*(8/100); % 0-1 --> oktas

% Convert longitudes to -180 to 180 range and sort to have monotonicaly
% increasing values
lon(lon > 180) = lon(lon > 180) - 360;
[lon_sort, sortIdx] = sort(lon);
cloudoktas_sort = cloudoktas(sortIdx,:,:);

% Swap lon and lat dimensions to get lat x lon x time
cloudoktas_sort_perm = permute(cloudoktas_sort, [2, 1, 3]);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CHECK VALUES AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(cloudoktas_sort_perm(:),100);

% Save the data
cloudcover = cloudoktas_sort_perm;
cloudcover_lat = lat;
cloudcover_lon = lon_sort;
save(fullpathOutputCloudFile,'cloudcover','cloudcover_lat','cloudcover_lon')

% Visual inspection
prepareDataForPlotting(fullpathOutputCloudFile,[],'oktas',...
    4,8,true,'fig_monthly_cloudcover_pincus','Cloud cover Pincus 2008')
