
% ======================================================================= %
%                                                                         %
%        Seawater dynamic viscosity climatology calculated from           %
%            salinity and temperature from WOA23 and the                  %
%              MIT seawater properties library routines                   %
%                                                                         %
% This script creates a global gridded climatology of seawater dynamic    %
% viscosity using salinity and temperature products from World Ocean      %
% Atlas 2023 and the MIT seawater properties library routines             %
% (https://web.mit.edu/seawater/). The output array is                    %
% 180 x 360 x 102 x 12 and has units of g cm-1 s-1.                       %                                        
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 10 Mar 2025                                   %
%                                                                         %
% ======================================================================= %

% Clear workspace, close figures, and add paths to plotting resources
close all; clear all; clc
addpath(genpath(fullfile('code')))
addpath(genpath(fullfile('resources','external','subaxis'))); 
addpath(genpath(fullfile('resources','external','subaxis')));
addpath(genpath(fullfile('resources','internal')));
addpath(genpath(fullfile('figures')))

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Output file
fullpathOutputVisco = fullfile('data','processed','dynamicvisco_calculated_woa23.mat');

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
% SECTION 3 - CALCULATE DYNAMIC VISCOSITY
% -------------------------------------------------------------------------

nLats = numel(woa_lat);  
nLons = numel(woa_lon);   
nDepths = numel(woa_depth_temp); 

% Units conversion of depth (from m to dbar)
pressure = NaN(nLats,nLons,nDepths);
for iLat = 1:nLats
    for iLon = 1:nLons
        pressure(iLat,iLon,:) = gsw_p_from_z(-1.*woa_depth_temp, woa_lat(iLat)); % Depth is negative in TEOS-10
    end
end

% Units conversion of salinity (from PSU to absolute units (g/kg))
salinityAbs = NaN(size(sal));
for iMonth = 1:12
    for iDepth = 1:size(sal,3)
        salinityAbs(:,:,iDepth,iMonth) = gsw_SA_from_SP(squeeze(sal(:,:,iDepth,iMonth)),...
            squeeze(pressure(:,:,iDepth)),woa_lon,woa_lat);
    end
end

% Calculation of dynamic viscosity 
dynvisco = NaN(size(temp));
for iMonth = 1:12
    for iDepth = 1:size(sal,3)
        dynvisco(:,:,iDepth,iMonth) = SW_Viscosity(temp(:,:,iDepth,iMonth),'C',...
            salinityAbs(:,:,iDepth,iMonth),'ppt'); % kg m-1 s-1
    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CHECK AND SAVE THE DATA
% -------------------------------------------------------------------------

% Check for spurious data points
% figure(); histogram(dynvisco(:), 100);

% Save output
dynvisco_lat = woa_lat;
dynvisco_lon = woa_lon;
dynvisco_depth = woa_depth_temp;
save(fullfile(fullpathOutputVisco),'dynvisco','dynvisco_lat','dynvisco_lon','dynvisco_depth','-v7.3')   

% Visual inspection
prepareDataForPlotting(fullpathOutputVisco,11,'kg m^{-1} s^{-1}',...
    5e-4,3e-3,true,'fig_monthly_dynvisco_calculated_woa23','Dynamic viscosity at 50 m, calculated from WOA23')
