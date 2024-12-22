
% ======================================================================= %
%                                                                         %
%                   kd climatology from NASA Aqua-MODIS                   % 
%                                                                         %
% This script reads in diffuse attenuation coefficient at 490 nm (kd)     %
% data product from NASA's Ocean Color website.                           %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://oceancolor.gsfc.nasa.gov/l3/order/                          % 
%   - Time resolution: L3 monthly climatology                             %
%   - Space resolution: 4 km (8640 x 4320)                                %
%   - Version: 2002–2024 (single sensor, Aqua-MODIS)                      %
%   - Units: m-1                                                          %
%                                                                         % 
% This script uses this external function:                                %
%       processSensorDataFromNASA - custom function                       %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 12 Nov 2024                                   %
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

fullpathInputKdAquaModisDir = fullfile('data','raw','kd','AquaMODIS');
fullpathOutputKdAquaModisFile = fullfile('data','processed','kd_aquamodis.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - READ IN KD FROM AQUA-MODIS
% -------------------------------------------------------------------------

% Get information
S = ncinfo(fullfile(fullpathInputKdAquaModisDir,... 
    ['AQUA_MODIS.' '20020701_20240731.L3m.MC.KD.Kd_490.4km' '.nc']));

% Read in the data
[kd,lat,lon] = processSensorDataFromNASA(fullpathInputKdAquaModisDir,...
    'AQUA_MODIS.','Kd_490','20020701_20240731.L3m.MC.KD.Kd_490.4km');

% Check for spurious data points
figure(); histogram(kd(:), 100);

% Save the data
kd_lat = lat;
kd_lon = lon;
save(fullpathOutputKdAquaModisFile,'kd','kd_lat','kd_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputKdAquaModisFile,[],'m^{-1}',...
    0,0.2,true,'fig_monthly_kd_aquamodis','kd Aqua-MODIS')
