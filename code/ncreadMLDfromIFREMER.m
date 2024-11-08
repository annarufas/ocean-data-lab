
% ======================================================================= %
%                                                                         %
%                      MLD climatology from IFREMER                       %
%                                                                         %
% This script reads in mixed layer depth (MLD) data from IFREMER, using   % 
% MLD defined by a density criterion corresponding to a temperature       %
% decrease of 0.2°C.                                                      %
%                                                                         %
% Dataset characteristics:                                                %
%   - http://www.ifremer.fr/cerweb/deboyer/mld/Surface_Mixed_Layer_Depth.php
%   - Time resolution: monthly climatology                                %
%   - Space resolution: 2° x 2° (180 x 90)                                %
%   - Version: dataset file was last updated on Nov 2008                  %
%   - Units: m                                                            %
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

filenameInputMld = 'mld_DR003_c1m_reg2.0.nc';
filenameOutputMld = 'mld_ifremer.mat';
fullpathInputMldFile = fullfile('data','raw','MLD','IFREMER',filenameInputMld);
fullpathOutputMldFile = fullfile('data','processed',filenameOutputMld);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

S = ncinfo(fullpathInputMldFile); % short summary

% Read in data
lat = double(ncread(fullpathInputMldFile,'lat'));
lon = double(ncread(fullpathInputMldFile,'lon'));
time = double(ncread(fullpathInputMldFile,'time'));
mixedlayerdepth = double(ncread(fullpathInputMldFile,'mld'));

% Some arrangements
mixedlayerdepth(mixedlayerdepth < 0)    = NaN; % -9999 represents missing data
mixedlayerdepth(mixedlayerdepth == 1e9) = NaN; % 1e9 represents land

% Sort longitudes to have monotonically increasing values
lonW = [(0:2:178),(-180:2:-2)]'; % combine positive and negative longitudes
[lonW_sort, sortIdx] = sort(lonW);
mld_sort = mixedlayerdepth(sortIdx,:,:);

% Swap lon and lat dimensions to get lat x lon x time
mld_sort_perm = permute(mld_sort, [2, 1, 3]);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - INSPECT AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(mld_sort_perm(:),100);

% Save the data
mld = mld_sort_perm;
mld_lat = lat;
mld_lon = lonW_sort;
save(fullpathOutputMldFile,'mld','mld_lat','mld_lon')

% Visual inspection
prepareDataForPlotting(fullpathOutputMldFile,[],'m',...
    0,200,true,'fig_monthly_mld_ifremer','MLD IFREMER')
