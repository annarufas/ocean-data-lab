
% ======================================================================= %
%                                                                         %
%          Zooplankton concentration climatology from CMIP6               %
%                                                                         %
% This script processes mesozooplankton concentration data (mol C m-3)    %
% from CMIP6, available through the Earth System Grid Federation (ESGF).  %
% (https://esgf-node.llnl.gov/projects/cmip6/). It reads .nc files        %
% previously downloaded from ESGF via 'wget', The standard CMIP variable  % 
% name for zooplankton concentration (zmeso) was identified at:           %
% https://clipc-services.ceda.ac.uk/dreq/mipVars.html. Models included in %
% this analysis: IPSL-PISCES, GFDL-COBALT and UKESM-MEDUSA.               %
%                                                                         %
% This script uses the external MATLAB function that is run in SLURM      %
% "regridZooplanktonConcentrationFromCMIP6.m". It regrids the data        %
% from a tripolar (native) grid to a regular latitude/longitude grid.     %
%                                                                         %
% Dataset characteristics:                                                %                                             
%   - http://esgf-node.llnl.gov/search/cmip6/?mip_era=CMIP6&activity_id=CMIP&frequency=mon&experiment_id=historical/ 
%   - Time resolution: monthly (1850-2014)                                %
%   - Space resolution: approx. 100 km                                    %
%   - Version: CMIP6                                                      %
%   - Units: mol C m-3 (input), mg C m-3 (output)                         % 
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

% Directory declarations
fullpathInputDataDir = fullfile('data','raw','zooplankton','CMIP6');
fullpathOutputDataDir = fullfile('data','processed');

% Input file declarations
filenameInputPisces  = 'zmeso_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_195001-201412.nc';
filenameInputCobalt = {'zmeso_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.nc',...
                       'zmeso_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_201001-201412.nc'};
filenameInputMedusa  = 'zmeso_Omon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-199912.nc';

% Ouput file declarations 
fullpathPiscesOutputFile = fullfile(fullpathOutputDataDir,'mesozoo_cmip6_pisces.mat');
fullpathCobaltOutputFile = fullfile(fullpathOutputDataDir,'mesozoo_cmip6_cobalt.mat');
fullpathMedusaOutputFile = fullfile(fullpathOutputDataDir,'mesozoo_cmip6_medusa.mat');
fullpathIntermediateFile = fullfile(fullpathInputDataDir,'cmip6_mesozooplankton_raw.mat');

% Grid that will be used to regrid data from high to lower resolution
fullpathGridFile = fullfile('data','raw','grid','grid_MITgcm_2p8deg.mat'); % MITgcm ocean model grid

% Set the target start date for data extraction
targetStartDates = struct('PISCES', datetime(1997,1,1), ...
                          'COBALT', datetime(1997,1,1), ...
                          'MEDUSA', datetime(1982,1,1));

% Parameters
MOLAR_MASS_CARBON = 12.011; % g C mol-1
NUM_DEPTH_ZONES = 5; % number of depth zones for regridding

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - PROGRAMMATICALLY DOWNLOAD .NC FILES CMIP6 USING 'WGET'
% -------------------------------------------------------------------------

% To download the .nc files, refer to "README_wget.txt" located in 
% ./data/raw/zooplankton/CMIP6/. This file provides instructions for 
% executing wget scripts to retrieve and store the required data files.

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - READ THE .NC FILES DOWNLOADED 
% -------------------------------------------------------------------------

% PISCES (uses native grid, 'gn')

fullpathPiscesFile = fullfile(fullpathInputDataDir,filenameInputPisces);
pisces.info  = ncinfo(fullpathPiscesFile); 
pisces.lat   = ncread(fullpathPiscesFile,'nav_lat');
pisces.lon   = ncread(fullpathPiscesFile,'nav_lon');
pisces.depth = ncread(fullpathPiscesFile,'olevel');

% Convert time into datetime units and find a start time
pisces_time_tmp = ncread(fullpathPiscesFile,'time');
pisces_time_datetime = datetime(1850,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(pisces_time_tmp);
[~,idxStartDate] = min(abs(pisces_time_datetime - targetStartDates.PISCES));
pisces.time = pisces_time_datetime(idxStartDate:end);

% Read only data for the time of interest
coordDims = size(pisces.lat);
startIndices = [1 1 1 idxStartDate]; 
dimCounts = [coordDims(1) coordDims(2) length(pisces.depth) length(pisces.time)]; 
pisces.mesozoo = ncread(fullpathPiscesFile,'zmeso',startIndices,dimCounts); % 362 x 332 x 75 x time

% .........................................................................

% COBALT (uses regular grid, 'gr')

fullpathCobaltFile = fullfile(fullpathInputDataDir,filenameInputCobalt);
cobalt.info  = ncinfo(fullpathCobaltFile{1}); 
cobalt.lat   = ncread(fullpathCobaltFile{1},'lat');
cobalt.lon   = ncread(fullpathCobaltFile{1},'lon'); % this needs to be rearranged as it goes from 0 to 360
cobalt.depth = ncread(fullpathCobaltFile{1},'lev');

% Convert time into datetime units and find a start time
for iFile = 1:numel(filenameInputCobalt)
    cobalt_time_tmp = ncread(fullpathCobaltFile{iFile},'time');
    cobalt_time_datetime = datetime(1850,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(cobalt_time_tmp);
    if (iFile == 1)
        [~,idxStartDate] = min(abs(cobalt_time_datetime - targetStartDates.COBALT));
        cobalt_time_first = cobalt_time_datetime(idxStartDate:end);
        startIndices = [1 1 1 idxStartDate]; 
        dimCounts = [numel(cobalt.lon) numel(cobalt.lat) numel(cobalt.depth) numel(cobalt_time_first)];
        cobalt_mesozoo_first = ncread(fullpathCobaltFile{iFile},'zmeso',startIndices,dimCounts);
    else
        cobalt_time_second = cobalt_time_datetime(:);
        startIndices = [1 1 1 1]; 
        dimCounts = [numel(cobalt.lon) numel(cobalt.lat) numel(cobalt.depth) numel(cobalt_time_second)];
        cobalt_mesozoo_second = ncread(fullpathCobaltFile{iFile},'zmeso',startIndices,dimCounts);
    end
end

% Combine time data
cobalt.time = [cobalt_time_first; cobalt_time_second];        
cobalt_mesozoo = cat(4, cobalt_mesozoo_first, cobalt_mesozoo_second); % 360 x 180 x 35 x time

% Adjust longitude from 0-360 to -180 to 180
cobalt.lon(cobalt.lon > 180) = cobalt.lon(cobalt.lon > 180) - 360; 

% Sort longitude to have monotonically increasing values
[lon_sorted, sortIdx] = sort(cobalt.lon);
cobalt.mesozoo = cobalt_mesozoo(sortIdx,:,:,:);
cobalt.lon = lon_sorted;

% .........................................................................

% MEDUSA (uses native grid, 'gn')

fullpathMedusaFile = fullfile(fullpathInputDataDir,filenameInputMedusa);
medusa.info  = ncinfo(fullpathMedusaFile); 
medusa.lat   = ncread(fullpathMedusaFile,'latitude');
medusa.lon   = ncread(fullpathMedusaFile,'longitude');
medusa.depth = ncread(fullpathMedusaFile,'lev');

% Convert time into datetime units and find a start time
medusa_time_tmp = ncread(fullpathMedusaFile,'time');
medusa_time_datetime = datetime(1850,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(medusa_time_tmp);
[~,idxStartDate] = min(abs(medusa_time_datetime - targetStartDates.MEDUSA));
medusa.time = medusa_time_datetime(idxStartDate:end);

% Read only data for the time of interest
coordDims = size(medusa.lat);
startIndices = [1 1 1 idxStartDate]; 
dimCounts = [coordDims(1) coordDims(2) length(medusa.depth) length(medusa.time)]; 
medusa.mesozoo = ncread(fullpathMedusaFile,'zmeso',startIndices,dimCounts); % 360 x 330 x 75 x time

% Realise that MEDUSA contains negative values, set those to 0
medusa.mesozoo(medusa.mesozoo < 0) = 0;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - SAVE THE DATA IN THEIR ORIGINAL GRIDS
% -------------------------------------------------------------------------

save(fullpathIntermediateFile,'pisces','cobalt','medusa','-v7.3')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - CREATE CLIMATOLOGIES AND PREPARE ARRAYS THAT NEED REGRIDDING
% -------------------------------------------------------------------------

% Create climatologies for each model. The PISCES and MEDUSA climatologies 
% will be regridded from native to regular grids in the next step. To 
% optimise regridding speed, these climatologies are divided into five 
% depth zones.

% PISCES
[pisces_climatology_native] = createClimatology(pisces,MOLAR_MASS_CARBON);
splitIntoDepthZones(pisces_climatology_native,pisces,NUM_DEPTH_ZONES,fullpathInputDataDir,'Pisces')

% COBALT
[mesozooClimatologyRegularCobalt] = createClimatology(cobalt,MOLAR_MASS_CARBON);
save(fullfile(fullpathInputDataDir,'mesozooClimatologyRegularCobalt.mat'),...
    'mesozooClimatologyRegularCobalt');

% MEDUSA
[medusa_climatology_native] = createClimatology(medusa,MOLAR_MASS_CARBON);
splitIntoDepthZones(medusa_climatology_native,medusa,NUM_DEPTH_ZONES,fullpathInputDataDir,'Medusa')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - RUN REGRIDDING FUNCTION IN SLURM
% -------------------------------------------------------------------------

% In SLURM, execute the script 'submit_zoo_regridding.sh', which will
% execute the MATLAB function 'regridZooplanktonConcentrationFromCMIP6.m'.
% Once finished, download the outputs to ./data/raw/zooplankton/CMIP6/'.

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - REPOSITION REGRIDDED DATA ON A LAT/LON ARRAY, SAVE AND
% VISUALISE ALL DATA
% -------------------------------------------------------------------------

saveClimatologiesInRegularGrid(fullpathIntermediateFile,...
    fullpathGridFile,fullpathInputDataDir,NUM_DEPTH_ZONES,...
    fullpathPiscesOutputFile,fullpathCobaltOutputFile,fullpathMedusaOutputFile)

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

% *************************************************************************

function [climatologyMesozooData] = createClimatology(mesozooStruct,MOLAR_MASS_CARBON)

    [nLon,nLat,nDepth,~] = size(mesozooStruct.mesozoo);

    % Loop over each month and calculate the climatology
    climatologyMesozooData = NaN(nLon,nLat,nDepth,12);
    months = month(mesozooStruct.time);
    for iMonth = 1:12
        idx = (months == iMonth);
        climatologyMesozooData(:,:,:,iMonth) =...
            mean(mesozooStruct.mesozoo(:,:,:,idx),4,'omitnan')*(MOLAR_MASS_CARBON*1e3); % mol C m-3 --> mg C m-3
    end

end % createClimatology

% *************************************************************************

function splitIntoDepthZones(mesozooClim,mesozooStruct,NUM_DEPTH_ZONES,...
    fullpathInputDataDir,modelName)

    [~,~,nDepths,~] = size(mesozooClim);
    lats = mesozooStruct.lat;
    lons = mesozooStruct.lons;

    depthZoneSize = floor(nDepths / NUM_DEPTH_ZONES);
    for iZone = 1:NUM_DEPTH_ZONES
        idxDepthStart = (iZone-1) * depthZoneSize + 1;
        if iZone == NUM_DEPTH_ZONES
            idxDepthEnd = nDepths;
        else
            idxDepthEnd = iZone * depthZoneSize;
        end
        mesozooClimatologyNative = mesozooClim(:,:,idxDepthStart:idxDepthEnd,:); % mg C m-3
        save(fullfile(fullpathInputDataDir,sprintf('mesozooClimatologyNative%s_dz%d.mat', modelName, iZone)),... 
            'mesozooClimatologyNative','lats','lons');
    end
    clear mesozooClimatologyNative lats lons; % clear temporary variables

end % splitIntoDepthZones

% *************************************************************************

function saveClimatologiesInRegularGrid(fullpathRawDataIntermediateFile,...
    fullpathGridFile,fullpathInputDataDir,NUM_DEPTH_ZONES,...
    fullpathPiscesOutputFile,fullpathCobaltOutputFile,fullpathMedusaOutputFile)

load(fullpathRawDataIntermediateFile,'pisces','cobalt','medusa')
load(fullpathGridFile,'Xbb','x','y','ixBb','iyBb'); 

nLocs = length(Xbb);
nRegriddedLat = length(y);
nRegriddedLon = length(x);

% .........................................................................

% PISCES

% Group depth zone data after regridding
mesozooClimatologyRegularPisces = NaN(nLocs,numel(pisces.depth),12); % 4448 x 75 x 12
depthZoneSize = floor(numel(pisces.depth) / NUM_DEPTH_ZONES);
for iZone = 1:NUM_DEPTH_ZONES
    load(fullfile(fullpathInputDataDir,sprintf('mesozooClimatologyRegularPisces_dz%d.mat',iZone)),...
        'mesozooClimatologyRegular');
    idxDepthStart = (iZone-1) * depthZoneSize + 1;
    if iZone == NUM_DEPTH_ZONES
        idxDepthEnd = numel(pisces.depth);
    else
        idxDepthEnd = iZone * depthZoneSize;
    end
    mesozooClimatologyRegularPisces(:,idxDepthStart:idxDepthEnd,:) = mesozooClimatologyRegular; % mg C m-3
    clear mesozooRegularPisces
end

% Reposition data on a lon x lat array
mesozoo_regrid = NaN(nRegriddedLon,nRegriddedLat,numel(pisces.depth),12);
for iLoc = 1:nLocs
    iLon = ixBb(iLoc);
    iLat = iyBb(iLoc);
    for iMonth = 1:12
        for iDepth = 1:numel(pisces.depth)
            mesozoo_regrid(iLon,iLat,iDepth,iMonth) = mesozooClimatologyRegularPisces(iLoc,iDepth,iMonth); % mg C m-3
        end
    end
end

% Swap lon and lat dimensions to get lat x lon x depth x time
mesozoo_regrid_perm = permute(mesozoo_regrid, [2, 1, 3, 4]);

% Save
mesozoo = mesozoo_regrid_perm;
mesozoo_lon = x;
mesozoo_lat = y;
mesozoo_depth = double(pisces.depth);
save(fullpathPiscesOutputFile,'mesozoo','mesozoo_lat','mesozoo_lon','mesozoo_depth');
clear mesozoo_lon mesozoo_lat mesozoo_depth mesozoo

% Visual inspection
prepareDataForPlotting(fullpathPiscesOutputFile,19,'mg C m^{-3}',...
    0,20,true,'fig_monthly_mesozoo_pisces','Mesozooplankton at 50 m, PISCES')

% .........................................................................

% COBALT

load(fullfile(fullpathInputDataDir,'mesozooClimatologyRegularCobalt.mat'),...
    'mesozooClimatologyRegularCobalt');

% Swap lon and lat dimensions to get lat x lon x depth x time
mesozoo_perm = permute(mesozooClimatologyRegularCobalt, [2, 1, 3, 4]);

% Save
mesozoo_lon = double(cobalt.lon);
mesozoo_lat = double(cobalt.lat);
mesozoo_depth = double(cobalt.depth);
mesozoo = mesozoo_perm;
save(fullpathCobaltOutputFile,'mesozoo','mesozoo_lat','mesozoo_lon','mesozoo_depth');
clear mesozoo_lon mesozoo_lat mesozoo_depth mesozoo

% Visual inspection
prepareDataForPlotting(fullpathCobaltOutputFile,5,'mg C m^{-3}',...
    0,20,true,'fig_monthly_mesozoo_cobalt','Mesozooplankton at 50 m, COBALT')

% .........................................................................

% MEDUSA

% Group depth zone data after regridding
mesozooClimatologyRegularMedusa = NaN(length(Xbb),numel(medusa.depth),12); % 4448 x 75 x 12
depthZoneSize = floor(numel(medusa.depth) / NUM_DEPTH_ZONES);
for iZone = 1:NUM_DEPTH_ZONES
    load(fullfile(fullpathInputDataDir,sprintf('mesozooClimatologyRegularMedusa_dz%d.mat',iZone)),...
        'mesozooClimatologyRegular');
    idxDepthStart = (iZone-1) * depthZoneSize + 1;
    if iZone == NUM_DEPTH_ZONES
        idxDepthEnd = numel(medusa.depth);
    else
        idxDepthEnd = iZone * depthZoneSize;
    end
    mesozooClimatologyRegularMedusa(:,idxDepthStart:idxDepthEnd,:) = mesozooClimatologyRegular; % mg C m-3
    clear mesozooRegularMedusa
end

% Reposition data on a lon x lat array
mesozoo_regrid = NaN(nRegriddedLon,nRegriddedLat,numel(medusa.depth),12);
for iLoc = 1:nLocs
    iLon = ixBb(iLoc);
    iLat = iyBb(iLoc);
    for iMonth = 1:12
        for iDepth = 1:numel(medusa.depth)
            mesozoo_regrid(iLon,iLat,iDepth,iMonth) = mesozooClimatologyRegularMedusa(iLoc,iDepth,iMonth); % mg C m-3
        end
    end
end

% Swap lon and lat dimensions to get lat x lon x depth x time
mesozoo_regrid_perm = permute(mesozoo_regrid, [2, 1, 3, 4]);

% Save
mesozoo = mesozoo_regrid_perm;
mesozoo_lon = x;
mesozoo_lat = y;
mesozoo_depth = double(medusa.depth);
save(fullpathMedusaOutputFile,'mesozoo','mesozoo_lat','mesozoo_lon','mesozoo_depth');
clear mesozoo_lon mesozoo_lat mesozoo_depth mesozoo

% Visual inspection
prepareDataForPlotting(fullpathMedusaOutputFile,19,'mg C m^{-3}',...
    0,20,true,'fig_monthly_mesozoo_medusa','Mesozooplankton at 50 m, MEDUSA')

end % saveClimatologiesInRegularGrid

% *************************************************************************