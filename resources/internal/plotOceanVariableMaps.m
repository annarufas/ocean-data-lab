function plotOceanVariableMaps(oceanVar,lonVector,latVector,myColormap,cbString,...
    caxisMin,caxisMax,isCommonColourBar,labelVars,figureName,sgString)

% PLOTOCEANVARIABLEMAPS Visualises oceanographic data across multiple 
% subplots, with the option for a common colour bar. It supports different 
% layouts based on the number of subplots.
%
%   INPUT: 
%       oceanVar          - variable to plot
%       lonVector         - vector of longitudes
%       latVector         - vector of latitudes
%       myColormap        - colour map of choice
%       cbString          - string label for the colour bar
%       caxisMin          - minimum value for colour axis scaling
%       caxisMax          - maximum value for colour axis scaling
%       isCommonColourBar - Boolean indicating whether to use a common colour bar for all subplots
%       figureName        - string specifying the base name for saved figure files
%       sgString          - string for a super title for the entire figure
%
% This script uses these external functions:
%       plotPcolorMap - custom function with instructions to plot
%       saveFigure    - custom function with type of figure output
%       m_map         - from FileExchange
%       subaxis       - from FileExchange
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

nSubplots = length(labelVars);

if (nSubplots == 4 && isCommonColourBar)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.45],'Color','w')
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(2,2,iSubplot,'Spacing',0.01,'Padding',0.02,'Margin',0.03);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColormap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',13)     

    end

    ax(3).pos(2) = ax(3).pos(2) + ax(1).pos(2)/6;
    ax(4).pos(2) = ax(4).pos(2) + ax(2).pos(2)/6;
    for iSubplot = 1:nSubplots
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 
    end

    cb = colorbar(haxis(4));
    cb.Location = 'southoutside';
    cb.Position(1) = 0.50-cb.Position(3)/2; 
    cb.Position(2) = ax(3).pos(2) - cb.Position(4);
    cb.Label.String = cbString;
    cb.FontSize = 12;

    saveFigure(figureName)

elseif (nSubplots == 4 && ~isCommonColourBar)
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.45],'Color','w')
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(2,2,iSubplot,'Spacing',0.01,'Padding',0.02,'Margin',0.04);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1 || iSubplot == 3)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.05;
        elseif (iSubplot == 2 || iSubplot == 4)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02;    
        end
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColormap,caxisMin(iSubplot),caxisMax(iSubplot));
        title(labelVars(iSubplot),'FontSize',13)  
        
        cb = colorbar('Location','eastoutside');
        cb.Position(1) = cb.Position(1) + 0.08;
        cb.Position(2) = cb.Position(2);
        cb.Position(3) = 0.015; % width
        cb.Position(4) = 0.22; % length
        cb.FontSize = 8;

    end

    saveFigure(figureName)
    
elseif (nSubplots == 3)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.50 0.28],'Color','w')
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(1,3,iSubplot,'Spacing',0.00,'Padding',0.005,'Margin',0.00);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.1;
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColormap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',13)     

    end

    cb = colorbar(haxis(2));
    cb.Location = 'southoutside';
    cb.Position(2) = cb.Position(2) - 0.10;  
    cb.Label.String = cbString;
    cb.FontSize = 12;

    saveFigure(figureName)

elseif (nSubplots == 5)
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.55],'Color','w')
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        if (iSubplot == 5) % combine positions 5 and 6 for bottom middle
            haxis(iSubplot) = subaxis(3,2,[5 6],'Spacing',0.01,'Padding',0.02,'Margin',0.04);
        else
            haxis(iSubplot) = subaxis(3,2,iSubplot,'Spacing',0.01,'Padding',0.02,'Margin',0.04);
        end 
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1 || iSubplot == 2)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.03;
        elseif (iSubplot == 3 || iSubplot == 4)    
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.045;
        elseif (iSubplot == 5)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.060;
        end
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColormap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',13)

    end
    
    cb = colorbar(haxis(5));
    cb.Location = 'southoutside';
    cb.Position(2) = cb.Position(2) - 0.09; %ax(3).pos(2) - cb.Position(4);
    cb.Label.String = cbString;
    cb.FontSize = 12;

    saveFigure(figureName)
 
elseif (nSubplots == 6 && ~isCommonColourBar)
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.55],'Color','w')
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(3,2,iSubplot,'Spacing',0.01,'Padding',0.02,'Margin',0.04);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1 || iSubplot == 3 || iSubplot == 5)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.05;
        elseif (iSubplot == 2 || iSubplot == 4 || iSubplot == 6)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02;    
        end
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColormap,caxisMin(iSubplot),caxisMax(iSubplot));
        title(labelVars(iSubplot),'FontSize',13)  
        
        cb = colorbar('Location','eastoutside');
        cb.Position(1) = cb.Position(1) + 0.08;
        cb.Position(2) = cb.Position(2);
        cb.Position(3) = 0.015; % width
        cb.Position(4) = 0.22; % length
        cb.FontSize = 8;

    end

    saveFigure(figureName)

elseif (nSubplots == 12)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.55 0.65],'Color','w')
    haxis = zeros(12,1);

    for iMonth = 1:12

        haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.005,'Padding',0.005,'Margin',0.085);
        ax(iMonth).pos = get(haxis(iMonth),'Position');
        ax(iMonth).pos(2) =  ax(iMonth).pos(2);
        set(haxis(iMonth),'Position',ax(iMonth).pos) 

        mydata = oceanVar(:,:,iMonth);
        plotPcolorMap(haxis(iMonth),lonVector,latVector,mydata,myColormap,caxisMin,caxisMax);
        title(labelVars(iMonth),'FontSize',12)     

    end
    
    sgtitle(sgString,'FontSize',14,'FontWeight','bold');

    cb = colorbar(haxis(12));
    cb.Location = 'southoutside';
    cb.Position(1) = 0.50-cb.Position(3)/2; 
    cb.Position(2) = ax(12).pos(2) - 0.02;
    cb.Label.String = cbString;
    cb.FontSize = 10;

    saveFigure(figureName)

end

end % plotOceanVariableMaps
