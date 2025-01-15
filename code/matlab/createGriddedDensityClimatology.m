
% ======================================================================= %
%                                                                         %
%               Seawater density climatology calculated from              %
%                   salinity and temperature from WOA23                   %
%                          and the SEAWATER toolbox                       %
%                                                                         %
% This script creates a global gridded climatology of seawater density    %
% using salinity and temperature products from World Ocean Atlas 2023     %
% and the SEAWATER toolbox. The output array is 180 x 360 x 102 x 12 and  % 
% has units of kg m-3.                                                    %                                        
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 8 Jan 2024                                    %
%                                                                         %
% ======================================================================= %

% Clear workspace, close figures, and add paths to plotting resources
close all; clear all; clc
addpath(genpath(fullfile('code')))
addpath(genpath(fullfile('resources','external'))); 
addpath(genpath(fullfile('resources','internal')));
addpath(genpath(fullfile('figures')))

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Output file
fullpathOutputRho = fullfile('data','processed','rho_calculated_woa23.mat');

% Input files
fullpathInputTempWoa = fullfile('data','processed','temp_monthly_woa23.mat');
fullpathInputSalWoa = fullfile('data','processed','sal_monthly_woa23.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD INPUT DATA
% -------------------------------------------------------------------------

load(fullpathInputTempWoa,'temp','woa_lat','woa_lon','woa_depth_temp');
load(fullpathInputSalWoa,'sal'); % same lat, lon and depths as temp

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE DENSITY
% -------------------------------------------------------------------------

% Prepare depth array
depth3 = repmat(woa_depth_temp,[1 length(woa_lat) length(woa_lon) 12]);
depth3 = permute(depth3,[2 3 1 4]);

% Calculate density using the SEAWATER toolbox
rho = sw_dens(sal,temp,depth3); % kg m-3

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CHECK AND SAVE THE DATA
% -------------------------------------------------------------------------

% Check for spurious data points
% figure(); histogram(rho(:), 100);

% Save output
rho_lat = woa_lat;
rho_lon = woa_lon;
rho_depth = woa_depth_temp;
save(fullfile(fullpathOutputRho),'rho','rho_lat','rho_lon','rho_depth','-v7.3')   

% Visual inspection
prepareDataForPlotting(fullpathOutputRho,11,'kg m^{-3}',...
    1020,1030,true,'fig_monthly_rho_calculated_woa23','Density at 50 m, calculated from WOA23')
