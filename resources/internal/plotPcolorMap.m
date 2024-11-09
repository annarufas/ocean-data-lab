function plotPcolorMap(haxis,lonVector,latVector,myData,myColormap,caxisMin,caxisMax)

    % Set projection
    m_proj('Robinson', 'clongitude', -63);

    % Swap latitude and longitude so that nRows = nLat and nCols = nLon
    [nRows,nCols] = size(myData);
    if (nRows > nCols)
        myData = permute(myData, [2, 1]);
    end

    myData(myData == Inf) = 0;

    % Because of the quirks of flat pcolor, there is a white line in the
    % join -180 to 180 (all the global projections have 360 deg ambiguities)
    m_pcolor(lonVector,latVector,myData);
    if any(lonVector(:) < 0)
        hold on;
        m_pcolor(lonVector-360,latVector,myData);
    end
    
    caxis([caxisMin caxisMax]);
%     shading flat;
    colormap(haxis,myColormap);

    m_grid('xtick', [-180:60:180], 'ytick', [-90:30:90], 'xticklabels', [], 'yticklabels', []); % , 'ticklen', 0.1, 'linestyle', 'none'
    m_coast('patch', [.7 .7 .7], 'edgecolor', 'black');
        
end % plotPcolorMap