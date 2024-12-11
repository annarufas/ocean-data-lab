
% ======================================================================= %
%                                                                         %
%            Aerosol dust deposition climatology from CMIP6               % 
%                                                                         %
% This script processes aerosol dust deposition data (kg m-2 s-1)         %
% downloaded from CMIP6 and available through the Earth System Grid       % 
% Federation (ESGF, https://esgf-node.llnl.gov/projects/cmip6/). It       % 
% reads in the .nc file previously downloaded from ESGF using a 'wget'    %
% function. The standard CMIP variable name for dust deposition (depdust) %
% was identified at: https://clipc-services.ceda.ac.uk/dreq/mipVars.html. %
% I have chosen NCAR-CESM2 historical run.                                %
%                                                                         %
% Dataset characteristics:                                                %                                             
%   - http://esgf-node.llnl.gov/search/cmip6/?mip_era=CMIP6&activity_id=CMIP&variable_id=depdust&frequency=mon&experiment_id=historical/ 
%   - Time resolution: monthly (1850-2014)                                %
%   - Space resolution: 100 km (288 x 192)                                %
%   - Version: CMIP6                                                      %
%   - Units: kg m-2 s-1 (input), g m-2 s-1 (output)                       % 
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

filenameInputDust = 'depdust_Emon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc';
fullpathInputDustFile = fullfile('data','raw','dust','CMIP6',filenameInputDust);
fullpathOutputDustFile = fullfile('data','processed','dustflux_cmip6_ncarcesm2.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - READ IN AEROSOL DUST DEPOSITION DATA
% -------------------------------------------------------------------------

S = ncinfo(fullpathInputDustFile); 

% Read in data
lat = double(ncread(fullpathInputDustFile,'lat'));
lon = double(ncread(fullpathInputDustFile,'lon'));
time = double(ncread(fullpathInputDustFile,'time'));
depdust = single(ncread(fullpathInputDustFile,'depdust')); % kg m-2 s-1

% Convert time into datetime units
refDate = datetime(1, 1, 1); % start date 'days since 0001-01-01 00:00:00' 
timeCalendar = refDate + days(time);

% Focus on the present-day (1985–2014) period
startDatePresent = datetime(1985, 1, 1);
endDatePresent = datetime(2014, 12, 31);
idxPresent = timeCalendar >= startDatePresent & timeCalendar <= endDatePresent;
timePresent = timeCalendar(idxPresent);
depdustPresent = depdust(:,:,idxPresent);

% Avoid negative flux numbers that will corrupt mean calculations
depdustPresent = abs(depdustPresent); 
       
% Create the climatology array
Dclim = NaN(numel(lon),numel(lat),12);
months = month(timePresent);
for iMonth = 1:12
    idx = (months == iMonth);
    Dclim(:,:,iMonth) = mean(depdustPresent(:,:,idx),3,'omitnan');
end

% Sort longitude
lon(lon > 180) = lon(lon > 180) - 360; % convert to -180 to 180
[lon_sort, sortIdx] = sort(lon);
Dclim_sort = Dclim(sortIdx,:,:);

% Swap lon and lat dimensions to get lat x lon x time
Dclim_sort_perm = permute(Dclim_sort, [2, 1, 3]);

% Units conversion
dustflux = Dclim_sort_perm.*1e3; % kg m-2 s-1 --> g m-2 s-1

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CHECK VALUES AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(dustflux(:),100);

% Save the data
dustflux_lon = lon_sort;
dustflux_lat = lat;
save(fullpathOutputDustFile,'dustflux','dustflux_lat','dustflux_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputDustFile,[],'g m^{-2} s^{-1}',...
    1e-10,1e-7,true,'fig_monthly_dust_cmip6','Aerosol dust deposition NCAR-CESM2')
 