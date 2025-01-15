
% ======================================================================= %
%                                                                         %
%             Carbonate ion, omega calcite and omega aragonite            %
%                     climatologies from GLODAP/CO2SYS                    %    
%                                                                         %
% This script processes global climatological data from the Global Ocean  %
% Data Analysis Project Version 2 (GLODAPv2.2016b) to compute carbonate   % 
% system variables. The calculations are performed using "CO2SYS.m", a    %
% program developed for CO2 system calculations (Lewis and Wallace, 1998; % 
% Van Heuven et al., 2011) and which you will find in the folder          %
% ./resources/external/.                                                  %
%                                                                         %
% Dataset characteristics:                                                %
%   - Data from: https://www.nodc.noaa.gov/archive/arc0107/0162565/2.2/data/0-data/mapped/                                  
%   - CO2SYS: https://www.nodc.noaa.gov/ocads/oceans/CO2SYS/co2rprt.html  %
%   - Time resolution: annual climatology                                 %
%   - Space resolution: 1° x 1° (360 x 180 x 33)                          %
%   - Version: v2, 2016                                                   %
%   - Units: Carbonate ion concentration: umol kg-1                       % 
%            Omega calcite: unitless                                      %
%            Omega aragonite: unitless                                    % 
%                                                                         %          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 7 Nov 2024                                    %
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

filenameInputCommonTag = 'GLODAPv2.2016b.';
fullpathInputDataDir = fullfile('data','raw','GLODAPv2');
fullpathOutputCarbonateIonDir = fullfile('data','processed','co3ion_co2sys.mat');
fullpathOutputOmegaCalciteDir = fullfile('data','processed','omegacalcite_co2sys.mat');
fullpathOutputOmegaAragoniteDir = fullfile('data','processed','omegaaragonite_co2sys.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

% These are the variables needed to run the formula that computes
% CO3= concentration, omega calcite and omega aragonite
varNames = {'TAlk', 'TCO2', 'salinity', 'silicate', 'PO4', 'temperature'};
nVars = numel(varNames);
data = cell(1, nVars);

% Load data from NetCDF files
for i = 1:nVars
    data{i} = ncread(fullfile(fullpathInputDataDir, [filenameInputCommonTag varNames{i} '.nc']), varNames{i});
end

% Extract variables
[talk, tco2, salinity, silicate, phosphate, temperature] = deal(data{:});
lon = ncread(fullfile(fullpathInputDataDir, [filenameInputCommonTag 'temperature.nc']), 'lon');
lat = ncread(fullfile(fullpathInputDataDir, [filenameInputCommonTag 'temperature.nc']), 'lat');
depth = ncread(fullfile(fullpathInputDataDir, [filenameInputCommonTag 'temperature.nc']), 'Depth');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CO2SYS CALCULATIONS
% -------------------------------------------------------------------------

% Initialise output arrays
co3 = NaN(numel(lon),numel(lat),numel(depth)); % umol kg-1
omegac = NaN(numel(lon),numel(lat),numel(depth)); % unitless
omegaa = NaN(numel(lon),numel(lat),numel(depth)); % unitless

% Parameter choices
par1type = 1; % The first parameter supplied is of type "1", which is "alkalinity"
par2type = 2; % The first parameter supplied is of type "2", which is "DIC"
pHscale  = 1; % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c    = 4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    = 1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

for iLon = 1:numel(lon)
    for iLat = 1:numel(lat)
        
        par1    = squeeze(talk(iLon,iLat,:)); 
        par2    = squeeze(tco2(iLon,iLat,:)); 
        sal     = squeeze(salinity(iLon,iLat,:)); 
        sil     = squeeze(silicate(iLon,iLat,:)); 
        po4     = squeeze(phosphate(iLon,iLat,:)); 
        tempin  = squeeze(temperature(iLon,iLat,:)); % Temperature at input conditions
        presin  = depth(:); % Pressure at input conditions
        tempout = squeeze(temperature(iLon,iLat,:)); % Temperature at output conditions - same as at input
        presout = depth(:); % Pressure at output conditions - same as at input
        
        [A,headers] = CO2SYS(par1,par2,par1type,par2type,sal,tempin,...
            tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
        co3(iLon,iLat,:) = A(:,22);
        omegac(iLon,iLat,:) = A(:,30);
        omegaa(iLon,iLat,:) = A(:,31);
        
    end
end

% Sort data
lonUnsorted = [(20.5:1:179.5), (-179.5:1:-0.5), (0.5:1:19.5)]'; % 20.5ºW to 19.5ºW
[lonSorted, idxSorting] = sort(lonUnsorted); % -179.5 to 179.5
co3Sorted = co3(idxSorting,:,:);
omegacSorted = omegac(idxSorting,:,:);
omegaaSorted = omegaa(idxSorting,:,:);

% Swap lon and lat dimensions to get lat x lon x depth
co3Perm    = permute(co3Sorted,    [2, 1, 3]);
omegacPerm = permute(omegacSorted, [2, 1, 3]);
omegaaPerm = permute(omegaaSorted, [2, 1, 3]);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - DATA VISUALISATION AND SAVING
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(co3Perm(:),10);
figure(); histogram(omegacPerm(:),10);
figure(); histogram(omegaaPerm(:),10);

% Save the data
glodap_lon = lonSorted;
glodap_lat = lat;
glodap_depth = depth;
carbonateion = co3Perm;
omegacalcite = omegacPerm;
omegaaragonite = omegaaPerm;
save(fullpathOutputCarbonateIonDir,'carbonateion','glodap_lat','glodap_lon','glodap_depth')
save(fullpathOutputOmegaCalciteDir,'omegacalcite','glodap_lat','glodap_lon','glodap_depth')
save(fullpathOutputOmegaAragoniteDir,'omegaaragonite','glodap_lat','glodap_lon','glodap_depth')

% Visual inspection
prepareDataForPlotting(fullpathOutputCarbonateIonDir,5,'umol kg^{-1}',...
    0,250,true,'fig_monthly_co3ion_co2sys','Carbonate ion at 50 m, CO2SYS')
prepareDataForPlotting(fullpathOutputOmegaCalciteDir,5,'unitless',...
    0,6,true,'fig_monthly_omegacalcite_co2sys','Omega calcite at 50 m, CO2SYS')
prepareDataForPlotting(fullpathOutputOmegaAragoniteDir,5,'unitless',...
    0,6,true,'fig_monthly_omegaaragonite_co2sys','Omega aragonite at 50 m, CO2SYS')
