function saveFigure(figureDirName,figureName)
    
    exportgraphics(gcf,fullfile('.','figures',figureDirName,...
        strcat(figureName,'.png')),'Resolution',600)
    
end