function prepareDataForPlotting(filePath,iDepth,cbString,caxisMin,caxisMax,...
    isCommonColourBar,figureName,sgString)

% PREPAREDATAFORPLOTTING Function designed to visualise monthly climatological
% oceanographic data from a specified file. If the size of the data exceeds
% a specified threshold set by MATLAB, the function regrids the data to 
% a lower resolution (1080 x 2160 pixels) using linear interpolation, so 
% that it is faster to plot. If the data are 4D (ie., include depth 
% dimension), the function extracts a specific depth slice based on the 
% input iDepth.
%
%   INPUT:
%       filePath          - path to the .mat file containing the data
%       iDepth            - index of the depth slice to visualise (if applicable)
%       cbString          - label for the colour bar
%       caxisMin          - minimum value for colour axis scaling
%       caxisMax          - maximum value for colour axis scaling
%       isCommonColourBar - flag to indicate whether to use a common colour bar for all subplots
%       figureName        - name of the figure file to be saved
%       sgString          - super title or general title for the figure
%
%   This script uses these external functions: 
%       brewermap             - from FileExchange
%       plotOceanVariableMaps - custom function
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 29 Oct 2024  
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

% Colour scale choice
myColourMap = flipud(brewermap(100,'RdYlBu')); % flip to have blue colours for low values

% Label subplots
labelMonths = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% Load the data
dataStruct = load(filePath);

% Get all variable names in the file
varNames = fieldnames(dataStruct);

% Identify lat, lon and data variables by name patterns
latVar = varNames{contains(varNames,'lat','IgnoreCase',true)};
lonVar = varNames{contains(varNames,'lon','IgnoreCase',true)};
depthIdx = contains(varNames,'depth','IgnoreCase',true);
if any(depthIdx)
    depthVar = varNames{depthIdx};
else
    depthVar = [];
end
errIdx = contains(varNames,'err','IgnoreCase',true);
if any(errIdx)
    errVar = varNames{errIdx};
else
    errVar = [];
end
dataVar = varNames{~contains(varNames,{'lat','lon','depth','err'},'IgnoreCase',true)};

% Extract lat, lon and data from the structure
lat = dataStruct.(latVar);
lon = dataStruct.(lonVar);
data = dataStruct.(dataVar);

% Regrid if the data array is too large and lon/lat come in the form of
% vectors (for already gridded lon/lat, don't apply this)
if min(size(lon)) == 1 && numel(lon) > 2160
    [lon, lat, data] = regridData(lon, lat, data, depthVar);
end

% Check if the data is 4D and, if so, select a depth slice
if ndims(data) == 4
    data = squeeze(data(:,:,iDepth,:)); % ith depth slice
end

% data(data == 0) = NaN; % to improve visualisation

% Check if the data spans multiple orders of magnitude (log10 values) and 
% decide whether to apply a log10 transformation. 
orderMin = log10(caxisMin);
if caxisMin == 0
    orderMin = 0;
end
orderMax = log10(caxisMax);
if (max(orderMax,orderMin) - min(orderMax,orderMin)) >= 4 % if the difference in orders of magnitude >= 4
    data = log10(data);
    caxisMin = log10(caxisMin);
    caxisMax = log10(caxisMax);
    cbString = [cbString, ' (log10)'];
end

% Plot
plotOceanVariableMaps(data,lon,lat,myColourMap,cbString,caxisMin,caxisMax,...
    isCommonColourBar,labelMonths,figureName,sgString)

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS TO THIS SCRIPT
% -------------------------------------------------------------------------

% *************************************************************************

function [newLon,newLat,regriddedData] = regridData(lon,lat,data,depthVar)

    newLon = linspace(-180, 180, 2160)';
    newLat = linspace(-90, 90, 1080)';

    % Determine lat/lon order before regridding
    nx = size(data, 1);
    ny = size(data, 2); 
    lonFirstDim = nx > ny; % 1 if longitude is the first dimension

    if isempty(depthVar) % if depth variable does not exist
        if lonFirstDim
            [X, Y, T] = ndgrid(lon, lat, (1:12)'); % indexes original grid
            [qX, qY, qT] = ndgrid(newLon, newLat, (1:12)'); % low-resolution grid
        else
            [X, Y, T] = ndgrid(lat, lon, (1:12)'); 
            [qX, qY, qT] = ndgrid(newLat, newLon, (1:12)'); 
        end
        F = griddedInterpolant(X, Y, T, data, 'linear', 'none');
        regriddedData = F(qX, qY, qT);
    elseif ~isempty(depthVar) % if depth variable exists
        if lonFirstDim
            [X, Y, Z, T] = ndgrid(lon, lat, (1:size(data, 3))', (1:12)'); 
            [qX, qY, qZ, qT] = ndgrid(newLon, newLat, (1:size(data, 3))', (1:12)'); 
        else
            [X, Y, Z, T] = ndgrid(lat, lon, (1:size(data, 3))', (1:12)'); 
            [qX, qY, qZ, qT] = ndgrid(newLat, newLon, (1:size(data, 3))', (1:12)'); 
        end
        F = griddedInterpolant(X, Y, Z, T, data, 'linear', 'none');
        regriddedData = F(qX, qY, qZ, qT);
    end

end % regridData

% *************************************************************************

end % prepareDataForPlotting
