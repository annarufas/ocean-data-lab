
% ======================================================================= %
%                                                                         %
%          NPP climatology from the Ocean Productivity Site               % 
%                                                                         %
% This script reads in net primary production (NPP) from Oregon State     %
% University's Ocean Productivity Site. The data are sourced from         %
% Aqua-MODIS and SeaWiFS sensors using Standard VGPM, CbPM2 and CAFE      %
% models.                                                                 %
%                                                                         %
% Dataset characteristics:                                                %
%   - http://orca.science.oregonstate.edu/npp_products.php                % 
%   - Time resolution: monthly climatology                                %
%   - Space resolution: 1/6 of a degree (1080 x 2160)                     %
%   - Version: 2024                                                       %
%   - Units: mg C m-2 d-1                                                 % 
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

fullpathInputCafeSeawifsDir = fullfile('data','raw','OceanProductivitySite','SeaWiFS','CAFE');
fullpathInputCafeModisDir = fullfile('data','raw','OceanProductivitySite','AquaMODIS','CAFE');
fullpathInputCbpmModisDir = fullfile('data','raw','OceanProductivitySite','AquaMODIS','CbPM');
fullpathInputVgpmModisDir = fullfile('data','raw','OceanProductivitySite','AquaMODIS','VGPM');

filenameOutputCafeSeawifs = 'npp_cafe_seawifs.mat';
filenameOutputCafeModis = 'npp_cafe_modis.mat';
filenameOutputCbpmModis = 'npp_cbpm_modis.mat';
filenameOutputVgpmModis = 'npp_vgpm_modis.mat';

fullpathOutputChlaFile = fullfile('data','processed');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - READ IN THE DATA
% -------------------------------------------------------------------------

% Process CAFE SeaWiFS
ncreadNPPdata(fullpathInputCafeSeawifsDir,'cafe.s',1998:2007,...
    filenameOutputCafeSeawifs,fullpathOutputChlaFile);

% Process CAFE MODIS 
ncreadNPPdata(fullpathInputCafeModisDir,'cafe.m',2003:2019,...
    filenameOutputCafeModis,fullpathOutputChlaFile);

% Process CbPM MODIS 
ncreadNPPdata(fullpathInputCbpmModisDir,'cbpm.m',2003:2019,...
    filenameOutputCbpmModis,fullpathOutputChlaFile);

% Process VGPM MODIS 
ncreadNPPdata(fullpathInputVgpmModisDir,'vgpm.m',2003:2019,...
    filenameOutputVgpmModis,fullpathOutputChlaFile);

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function ncreadNPPdata(inputDir, tagInputFilename, yearSequence, outputFilename, outputDir)
   
    tagParts = split(tagInputFilename, '.'); % splits at '.'

    % Extract HDF files from gz archives if needed
    hdfFiles = dir(fullfile(inputDir, '*.hdf'));
    if isempty(hdfFiles)
        for iYear = yearSequence
            untar(fullfile(inputDir, [tagInputFilename '.' sprintf('%04d',iYear) '.tar']));
        end
        gunzip('*.gz');
    end
 
    % Get the number of years
    nyrs = length(yearSequence);
    
    % Initialise the data array
    dataNpp = zeros(1080, 2160, 12, nyrs);
    
    % Loop through each year and month to process the data
    for iYear = 1:nyrs
        iTime = yearSequence(iYear);
        dd = dir(fullfile(inputDir, [tagParts{1} '.' num2str(iTime) '*.hdf']));
        days = zeros(12,1);
        for iMonth = 1:12
            ddtmp = split(dd(iMonth).name,'.');
            days(iMonth) = str2num(ddtmp{2});
        end
        [days,idx] = sort(days);
        dd = dd(idx);
        for iMonth = 1:12
             dataNpp(:,:,iMonth,iYear) = double(hdfread(fullfile(inputDir, dd(iMonth).name),'npp'));  
        end
    end
    
    % Replace negative values with NaN and compute the average
    dataNpp(dataNpp < 0) = NaN;
    dataNppAvg = mean(dataNpp, 4, 'omitnan');

    % Bin the data in months and calculate the standard deviation
    dataNppErr = NaN(1080,2160,12); % standard deviation
    dataNppN   = NaN(1080,2160,12);

    for iLat = 1:1080
        for iLon = 1:2160
            for iMonth = 1:12
                if (sum(~isnan(dataNpp(iLat,iLon,iMonth,:))) > 1)
                    dataNppErr(iLat,iLon,iMonth) = std(dataNpp(iLat,iLon,iMonth,:),'omitnan');
                    dataNppN(iLat,iLon,iMonth) = sum(dataNpp(iLat,iLon,iMonth,:) >= 0,'omitnan');
                elseif (sum(~isnan(dataNpp(iLat,iLon,iMonth,:))) == 1)
                    dataNppErr(iLat,iLon,iMonth) = 0;
                    dataNppN(iLat,iLon,iMonth) = 1;
                end
            end 
        end
    end

    % This is the lat and lon arrangement of the data
    lons = linspace(-180, 180, 2160);
    lats = linspace(90, -90, 1080);

    % Sort lat to have monotonically increasing values and apply the sorting to
    % the data
    [lats_sort, sortIdx] = sort(lats);
    npp_sort   = dataNppAvg(sortIdx,:,:);
    err_sort   = dataNppErr(sortIdx,:,:,:);

    % Save the data
    npp_avg = npp_sort;
    npp_err = err_sort;
    npp_lon = lons;
    npp_lat = lats_sort;
    save(fullfile(outputDir,outputFilename),'npp_avg','npp_err','npp_lon','npp_lat')
    
    % Visual inspection
    parts = strsplit(outputFilename, {'_', '.'});
    titlePart1 = upper(parts{2});  % 'vgpm' -> 'VGPM'
    titlePart2 = upper(parts{3});  % 'modis' -> 'MODIS'
    sgTitle = sprintf('%s %s %s', 'NPP', titlePart1, titlePart2);
    prepareDataForPlotting(fullfile(outputDir,outputFilename),[],'mg C m-2 d-1',...
        0,1000,true,strcat('fig_monthly_',strcat('npp_',parts{2},'_',parts{3})),sgTitle)

end % function ncreadNPPdata
