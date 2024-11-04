function plotPcolorMap(haxis,lonVector,latVector,myData,myColormap,caxisMin,caxisMax)

    m_proj('Robinson', 'clongitude', -63);
    
    [nRows,nCols] = size(myData);
    if (nRows > nCols)
        myData = flipud(rot90(myData));
    end
    
    myData(myData == Inf) = 0;
%     max(max(myData))
%     minValue = min(myData(myData > 0))
    
    m_pcolor(lonVector,latVector,myData);
    if any(lonVector(:) < 0)
        hold on;
        m_pcolor(lonVector-360,latVector,myData);
    end
    
    caxis([caxisMin caxisMax]);
    shading flat;
    colormap(haxis,myColormap);
    m_grid('xticklabels', [], 'yticklabels', [], 'ticklen', 0, 'linestyle', 'none');
    m_coast('patch', [.7 .7 .7], 'edgecolor', 'black');
        
end % plotPcolorMap