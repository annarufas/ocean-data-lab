function plotPcolorMap(haxis,lon,lat,myData,myColourMap,caxisMin,caxisMax)

% PLOTPCOLORMAP Function with my favourite settings to plot pcolor maps.
%
%   INPUT: 
%       haxis       - figure axis handle
%       lon         - longitudes (can be -180 to 180 or 0 to 360, vector or array)
%       lat         - latitudes (-90 to 90, vector or array)
%       myData      - 2D data array (can be lat x lon or lon x lat)
%       myColourMap - colour map of choice
%       caxisMin    - minimum value in the colour bar
%       caxisMax    - maximum value in the colour bar
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 9 Nov 2024 
%   Version 2.0 - 13 Jan 2025 (handles lat/lon in both vector and gridded form) 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

% Swap latitude and longitude so that nRows = nLat and nCols = nLon
[nRows,nCols] = size(myData);
if (min(size(lon)) == 1) && (nRows > nCols)
    myData = permute(myData, [2, 1]);
end

% Transform into gridded form if not already
if min(size(lon)) == 1 && min(size(lat)) == 1
    [X,Y] = meshgrid(lon,lat);
else
    X = lon;
    Y = lat;
end

% Set projection
m_proj('Robinson', 'clongitude', -63);

% Because of the quirks of flat pcolor, there is a white line in the
% join -180 to 180 (all the global projections have 360 deg ambiguities)
m_pcolor(X,Y,myData);
if any(X(:) < 0)
    hold on;
    m_pcolor(X-360,Y,myData);
end

% Set axis limits, colour scaling and grid
caxis([caxisMin caxisMax]);
shading flat;
colormap(haxis,myColourMap);
m_grid('xtick', [-180:60:180], 'ytick', [-90:30:90], 'xticklabels', [], 'yticklabels', []); % , 'ticklen', 0.1, 'linestyle', 'none'
m_coast('patch', [.7 .7 .7], 'edgecolor', 'black');

end % plotPcolorMap