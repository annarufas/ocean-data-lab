
% ======================================================================= %
%                                                                         %
%        PAR0 climatology calculated from trigonometric/astronomic        %
%         equations and input data of cloud cover and ice fraction        % 
%                                                                         %
% This script creates a global gridded climatology of photosynthetic      %
% available radiation at the surface ocean (PAR0) using trigonometric/    %
% astronomic equations, incorporating inputs of sea ice fraction (from    % 
% CMEMS) and cloud cover (from Pincus et al. (2008)). The output arrays   %
% are monthly (180 x 360 x 12 pixels) and daily (180 x 360 x 365 pixels)  %
% climatologies in units of W m-2.                                        %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
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

% Output files
fullpathOutputPar0dailyFile = fullfile('data','processed','par0_daily_calculated.mat');
fullpathOutputNumDaylightHoursFile = fullfile('data','processed','ndaylighthours_daily_calculated.mat');
fullpathOutputPar0monthlyFile = fullfile('data','processed','par0_monthly_calculated.mat');

% Cloud cover and ice fraction data
fullpathCloudCoverClimatology = fullfile('data','processed','cloudcover_pincus.mat');
fullpathIceFractionClimatology = fullfile('data','processed','icefrac_cmems_phys.mat');

% Grid used to extract ocean locations to calculate PAR0 (I don't want to 
% calculate PAR0 at each grid point in the lat/lon array, which includes 
% land locations)
fullpathGridFile = fullfile('data','raw','grid','grid_GEBCO_360_180.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD AUXILLIARY DATA
% -------------------------------------------------------------------------

% Load cloud cover and ice fraction data
load(fullpathCloudCoverClimatology,'cloudcover','cloudcover_lat','cloudcover_lon') % oktas
load(fullpathIceFractionClimatology,'icefrac','icefrac_lat','icefrac_lon') % 0-1 fraction

% Load the grid 
load(fullpathGridFile,'Xbb','Ybb','x','y','ixBb','iyBb') 
nOceanLocs = length(Ybb);
oceanLats = Xbb;
oceanLons = Ybb;
nOceanLats = length(x);
nOceanLons = length(y);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - INTERPOLATION FUNCTIONS
% -------------------------------------------------------------------------

% Grids for interpolation
[Xclo, Yclo, Tclo] = ndgrid(cloudcover_lat, cloudcover_lon, (1:12)');
[Xice, Yice, Tice] = ndgrid(icefrac_lat, icefrac_lon, (1:12)');

% Interpolation functions
Fclo = griddedInterpolant(Xclo, Yclo, Tclo, cloudcover, 'linear', 'none');
Fice = griddedInterpolant(Xice, Yice, Tice, icefrac, 'linear', 'none');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - DAILY PAR0 CALCULATIONS
% -------------------------------------------------------------------------

Pdaily = NaN(nOceanLocs,365);
nDaylightHours = NaN(nOceanLocs,365);

for iLoc = 1:nOceanLocs

    qLat = oceanLats(iLoc);
    qLon = oceanLons(iLoc);   

    avgPar0daylight = zeros(365,1); % W m-2
    totPar0daylight = zeros(365,1); % J m-2

    % Query points for interpolation
    [qX, qY, qT] = ndgrid(qLat, qLon, (1:365)');

    % Get cloud cover
    qCloudFrac = Fclo(qX, qY, qT);
    qCloudFrac(qCloudFrac<0) = 0;
    qCloudFrac(isnan(qCloudFrac)) = 0; 

    % Get ice fraction
    qIceFrac = Fice(qX, qY, qT);
    qIceFrac(qIceFrac<0) = min(qIceFrac(qIceFrac > 0)); % 0-1
    qIceFrac(isnan(qIceFrac)) = 0;

    % Calculate the average PAR0 received during the daylight period (nDaylightHours) of each day
    [avgPar0daylight, nDaylightHours(iLoc,:)] = calculatePAR0fromTrigonometricEquations(...
        qLat, squeeze(qCloudFrac), squeeze(qIceFrac)); % W m-2

    % Calculate total PAR0 during daylight and daily mean PAR0
    totPar0daylight = avgPar0daylight.*squeeze(nDaylightHours(iLoc,:))'*3600; % J m-2
    Pdaily(iLoc,:) = totPar0daylight./(24*3600); % J m-2 --> W m-2 (= J s-1 m-2)
    
end    

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - MONTHLY CLIMATOLOGY CALCULATIONS
% -------------------------------------------------------------------------

Pclim = NaN(nOceanLocs,12);
daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % days in each month, aligning with a non-leap year
idxDay = 1; 
for iMonth = 1:12
    currentMonthData = Pdaily(:,idxDay:idxDay+daysInMonth(iMonth)-1);
    Pclim(:,iMonth) = mean(currentMonthData, 2, 'omitnan');
    idxDay = idxDay + daysInMonth(iMonth); % move to the starting day of the next month
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - REPOSITION DATA ON A LON/LAT ARRAY
% -------------------------------------------------------------------------

geo_nDaylightHours = NaN(nOceanLats,nOceanLons,365);
geo_Pdaily = NaN(nOceanLats,nOceanLons,365);
geo_Pclim  = NaN(nOceanLats,nOceanLons,12);

for iLoc = 1:nOceanLocs
    iLon = iyBb(iLoc);
    iLat = ixBb(iLoc);
    geo_nDaylightHours(iLat,iLon,:) = nDaylightHours(iLoc,:);
    geo_Pdaily(iLat,iLon,:) = Pdaily(iLoc,:);
    geo_Pclim(iLat,iLon,:)  = Pclim(iLoc,:);
end 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - CHECKINGS AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
% figure(); histogram(geo_Pclim(:),100);
% figure(); histogram(geo_Pdaily(:),100);

% Save output
par0_lat = x';
par0_lon = y';
par0daily = geo_Pdaily;
par0 = geo_Pclim;
par0daylighthours = geo_nDaylightHours;
save(fullfile(fullpathOutputPar0dailyFile),'par0daily','par0_lat','par0_lon') 
save(fullfile(fullpathOutputNumDaylightHoursFile),'par0daylighthours','par0_lat','par0_lon') 
save(fullfile(fullpathOutputPar0monthlyFile),'par0','par0_lat','par0_lon') 

% Visual inspection
prepareDataForPlotting(fullpathOutputPar0monthlyFile,[],'W m^{-2}',...
    0,200,true,'fig_monthly_par0_calculated','PAR0 calculated from equations')
