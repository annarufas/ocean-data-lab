function saveFigure(figureName)
    
    exportgraphics(gcf,fullfile('.','figures',strcat(figureName,'.pdf')),...
        'Resolution',600)
    exportgraphics(gcf,fullfile('.','figures',strcat(figureName,'.png')),...
        'Resolution',600)  

end