function saveFigure(figureName)
    
    exportgraphics(gcf,fullfile('.','figures',strcat(figureName,'.pdf')),...
        'Resolution',600)
    
end