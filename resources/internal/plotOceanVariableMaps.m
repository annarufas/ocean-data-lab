function plotOceanVariableMaps(oceanVar,lonVector,latVector,myColourMap,cbString,...
    caxisMin,caxisMax,isCommonColourBar,labelVars,figureName,sgString)

% PLOTOCEANVARIABLEMAPS Visualises oceanographic data across multiple 
% subplots, with the option for a common colour bar. It supports different 
% layouts based on the number of subplots.
%
%   INPUT: 
%       oceanVar          - variable to plot
%       lonVector         - vector of longitudes
%       latVector         - vector of latitudes
%       myColourMap       - colour map of choice
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

% oceanVar = auxArray;
% lonVector = config.mappingProps.mapLons;
% latVector = config.mappingProps.mapLats;
% myColourMap = config.mappingProps.myColourMapOceanVars;
% isCommonColourBar = [];
% labelVars = titleStr;
% figureName = strcat('output_global_annual_sms_',auxName);
% sgString = [];
                     
% nSubplots = numel(labelVars);
nSubplots = size(oceanVar,3);

if (nSubplots == 1)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.35],'Color','w')
    
    [oceanVar,caxisMin,caxisMax,ticks,tickLabels] = assessLogarithmicTransformation(oceanVar);

    haxis = subaxis(1,1,1,'Spacing',0.01,'Padding',0.02,'Margin',0.03);
    ax = get(haxis,'Position'); % Store position for later adjustments

    mydata = oceanVar(:,:);
    plotPcolorMap(haxis,lonVector,latVector,mydata,myColourMap,caxisMin,caxisMax);
    title(labelVars,'FontSize',16)     

    cb = colorbar(haxis,'southoutside');
    cb.Label.String = cbString;
    
    % Set axis limits and ticks
    caxis([caxisMin caxisMax]);
    cb.Ticks = ticks;
    cb.TickLabels = tickLabels;
    cb.FontSize = 9; % font size of tick labels
    cb.Label.FontSize = 12; % font size of colour bar string

    saveFigure(figureName)
   
elseif (nSubplots == 4 && isCommonColourBar)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.45],'Color','w')
    haxis = zeros(nSubplots,1);
    
    if (isempty(caxisMax) && isempty(caxisMin))
        disp('isempty')
        [oceanVar,caxisMin,caxisMax,ticks,tickLabels] = assessLogarithmicTransformation(oceanVar);
    else
        disp('is not empty')
        % Determine step size dynamically (max 8 ticks)
        stepSizes = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
        stepSize = stepSizes(find((caxisMax - caxisMin) ./ stepSizes <= 8, 1, 'first'));

        % Generate tick positions with at least 3 ticks but no more
        % than 8
        ticks = caxisMin:stepSize:caxisMax;
        if numel(ticks) < 3
            ticks = linspace(caxisMin, caxisMax, 3);
        elseif numel(ticks) > 8
            ticks = linspace(caxisMin, caxisMax, 8);
        end
        
        % Format tick labels
        % Ensure caxisMin is not zero (log10(0) is undefined)
        if caxisMin == 0
            caxisMin = min(ticks(ticks > 0)); % Set to the smallest positive tick
        end
        if abs(log10(caxisMin)) >= 4 || abs(log10(caxisMax)) >= 4
            tickLabels = arrayfun(@(x) sprintf('%.1e', x), ticks, 'UniformOutput', false);
        else
            tickLabels = arrayfun(@(x) sprintf('%g', x), ticks, 'UniformOutput', false);
        end

    end

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(2,2,iSubplot,'Spacing',0.01,'Padding',0.02,'Margin',0.03);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.02; 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColourMap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',13)     

    end

    ax(3).pos(2) = ax(3).pos(2) + ax(1).pos(2)/5;
    ax(4).pos(2) = ax(4).pos(2) + ax(2).pos(2)/5;
    for iSubplot = 1:nSubplots
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 
    end

    cb = colorbar(haxis(4));
    cb.Location = 'southoutside';
    cb.Position(1) = 0.50-cb.Position(3)/2; 
    cb.Position(2) = ax(3).pos(2) - cb.Position(4)/2;
    cb.Label.String = cbString;
    
    % Set axis limits and ticks
    caxis([caxisMin caxisMax]);
    cb.Ticks = ticks;
    cb.TickLabels = tickLabels;
    cb.FontSize = 9; % font size of tick labels
    cb.Label.FontSize = 12; % font size of colour bar string

    % Give common title to your figure
    if ~isempty(sgString)
        sgtitle(sgString,'FontSize',16,'FontWeight','bold');
    end

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
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColourMap,caxisMin(iSubplot),caxisMax(iSubplot));
        title(labelVars(iSubplot),'FontSize',13)  
        
        cb = colorbar('Location','eastoutside');
        cb.Position(1) = cb.Position(1) + 0.08;
        cb.Position(2) = cb.Position(2);
        cb.Position(3) = 0.015; % width
        cb.Position(4) = 0.22; % length
        cb.FontSize = 7; % font size of tick labels
        cb.Label.FontSize = 8; % font size of colour bar string

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
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColourMap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',13)     

    end

    cb = colorbar(haxis(2));
    cb.Location = 'southoutside';
    cb.Position(2) = cb.Position(2) - 0.10;  
    cb.Label.String = cbString;
    cb.FontSize = 10; % font size of tick labels
    cb.Label.FontSize = 12; % font size of colour bar string

    saveFigure(figureName)

elseif (nSubplots == 5 && isCommonColourBar)
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.55],'Color','w')
    haxis = zeros(nSubplots,1);

    [oceanVar,caxisMin,caxisMax,ticks,tickLabels] = assessLogarithmicTransformation(oceanVar);
    
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
    
    % Set axis limits and ticks
    caxis([caxisMin caxisMax]);
    cb.Ticks = ticks;
    cb.TickLabels = tickLabels;
    cb.FontSize = 10; % font size of tick labels
    cb.Label.FontSize = 12; % font size of colour bar string

    saveFigure(figureName)

elseif (nSubplots == 6 && isCommonColourBar)
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.55],'Color','w')
    haxis = zeros(nSubplots,1);

    [oceanVar,caxisMin,caxisMax,ticks,tickLabels] = assessLogarithmicTransformation(oceanVar);
    
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
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColormap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',13)

    end
    
    cb = colorbar(haxis(6));
    cb.Location = 'southoutside';
    cb.Position(1) = 0.25; 
    cb.Position(3) = 0.50; % length 
    cb.Position(4) = 0.02; % width 
    cb.Label.String = cbString;
    
    % Set axis limits and ticks
    caxis([caxisMin caxisMax]);
    cb.Ticks = ticks;
    cb.TickLabels = tickLabels;
    cb.FontSize = 10; % font size of tick labels
    cb.Label.FontSize = 12; % font size of colour bar string

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
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColourMap,caxisMin(iSubplot),caxisMax(iSubplot));
        title(labelVars(iSubplot),'FontSize',13)  
        
        cb = colorbar('Location','eastoutside');
        cb.Position(1) = cb.Position(1) + 0.08;
        cb.Position(2) = cb.Position(2);
        cb.Position(3) = 0.015; % width
        cb.Position(4) = 0.22; % length
        cb.FontSize = 8;

    end

    saveFigure(figureName)
    
elseif (nSubplots == 8 && ~isCommonColourBar)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.55],'Color','w')
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(4,2,iSubplot,'Spacing',0.01,'Padding',0.01,'Margin',0.045);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1 || iSubplot == 3 || iSubplot == 5 || iSubplot == 7)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02;
        elseif (iSubplot == 2 || iSubplot == 4 || iSubplot == 6 || iSubplot == 8)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.01;    
        end
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.03;    
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos)

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColourMap,caxisMin(iSubplot),caxisMax(iSubplot));
        if (iSubplot == 1 || iSubplot == 2)
            tt = title(labelVars(iSubplot),'FontSize',12);
            tt.Position(2) = tt.Position(2) + 0.10;
        end
        
        % Annotation to the left of the subplot
        if (iSubplot == 1 || iSubplot == 3 || iSubplot == 5 || iSubplot == 7)
            
            % Define the position for the annotation relative to the subplot
            x_text = ax(iSubplot).pos(1) - 0.12; % Move left outside the subplot by 0.05 normalized units
            y_text = 0.50;% + ax(iSubplot).pos(4)/2; % Vertical center
    
            % Add vertical text to the left (westoutside)
            text(x_text, y_text, cbString(iSubplot), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', ...
                'Rotation', 90, 'Units', 'normalized');
    
        end
        
        cb = colorbar('Location','eastoutside');
        cb.Position(1) = cb.Position(1) + 0.06;
        cb.Position(2) = cb.Position(2) + 0.01;
        cb.Position(3) = 0.015; % width
        cb.Position(4) = 0.15; % length
        cb.FontSize = 8;
        
        if caxisMax(iSubplot) < 0 % negative number indicates log10 (e.g., log(1e-16)= -16)
            nTicks = (caxisMax(iSubplot) - caxisMin(iSubplot))+1;
            ticks = linspace(caxisMin(iSubplot),caxisMax(iSubplot),nTicks);
            % Adjusts the format dynamically, ensuring the number of decimals 
            % is appropriate based on the scale
            tickLabels = arrayfun(@(x)... 
                sprintf(['%.' num2str(max(0, -floor(log10(abs(10^x)))) ) 'f'], 10^x), ticks, 'UniformOutput', false);
            cb.Ticks = ticks; % set the ticks in original scale 
            cb.TickLabels = tickLabels; % set the custom tick labels
        end

    end
    
    % Give common title to your figure
    if ~isempty(sgString)
        sgtitle(sgString,'FontSize',14,'FontWeight','bold');
    end

    saveFigure(figureName)

elseif (nSubplots == 9 && isCommonColourBar)
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.55],'Color','w')
    haxis = zeros(nSubplots,1);

    [oceanVar,caxisMin,caxisMax,ticks,tickLabels] = assessLogarithmicTransformation(oceanVar);
    
    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(3,3,iSubplot,'Spacing',0.01,'Padding',0.02,'Margin',0.04);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1 || iSubplot == 2 || iSubplot == 3)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.02;
        elseif (iSubplot == 4 || iSubplot == 5 || iSubplot == 6)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.04;
        elseif (iSubplot == 7 || iSubplot == 8 || iSubplot == 9)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.06;
        end  
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColormap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',13)

    end
    
    cb = colorbar(haxis(nSubplots));
    cb.Location = 'southoutside';
    cb.Position(1) = 0.25; 
    cb.Position(2) = 0.075; 
    cb.Position(3) = 0.50; % length 
    cb.Position(4) = 0.02; % width 
    cb.Label.String = cbString;
    
    % Set axis limits and ticks
    caxis([caxisMin caxisMax]);
    cb.Ticks = ticks;
    cb.TickLabels = tickLabels;
    cb.FontSize = 9; % font size of tick labels
    cb.Label.FontSize = 11; % font size of colour bar string

    saveFigure(figureName)
    
elseif (nSubplots == 12)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.55 0.65],'Color','w')
    haxis = zeros(nSubplots,1);

    for iMonth = 1:12

        haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.005,'Padding',0.005,'Margin',0.085);
        ax(iMonth).pos = get(haxis(iMonth),'Position');
        ax(iMonth).pos(2) = ax(iMonth).pos(2) + 0.005;
        set(haxis(iMonth),'Position',ax(iMonth).pos) 

        mydata = oceanVar(:,:,iMonth);
        plotPcolorMap(haxis(iMonth),lonVector,latVector,mydata,myColourMap,caxisMin,caxisMax);
        title(labelVars(iMonth),'FontSize',12)     

    end
    
    sgtitle(sgString,'FontSize',14,'FontWeight','bold');
    
    cb = colorbar(haxis(12));
    cb.Location = 'southoutside';
    cb.Position(1) = 0.50-cb.Position(3)/2; 
    cb.Position(2) = ax(12).pos(2) - 0.02;
    
    if contains(cbString,'log10')
        nTicks = (caxisMax - caxisMin)+1;
        ticks = linspace(caxisMin,caxisMax,nTicks);
        % Adjusts the format dynamically, ensuring the number of decimals 
        % is appropriate based on the scale
        tickLabels = arrayfun(@(x)... 
            sprintf(['%.' num2str(max(0, -floor(log10(abs(10^x)))) ) 'f'], 10^x), ticks, 'UniformOutput', false);
        cb.Ticks = ticks; % set the ticks in original scale 
        cb.TickLabels = tickLabels; % set the custom tick labels
    end
    
    cbStringModified = strrep(cbString, '(log10)', ''); % remove 'log10' from the string
    cb.Label.String = cbStringModified; % set colorbar label without 'log10'
    cb.Label.Rotation = 0;
    cb.FontSize = 10;

    saveFigure(figureName)

elseif (nSubplots == 14 && ~isCommonColourBar)

    figure()
    set(gcf,'Units','Normalized','Position',[0.0 0.0 0.50 0.89],'Color','w') 
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(5,3,iSubplot,'Spacing',0.015,'Padding',0.015,'Margin',0.027);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1 || iSubplot == 4 || iSubplot == 7 || iSubplot == 10)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.03; 
        elseif (iSubplot == 2 || iSubplot == 5 || iSubplot == 8 || iSubplot == 11)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02; 
        elseif (iSubplot == 3 || iSubplot == 6 || iSubplot == 9 || iSubplot == 12)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.01; 
        elseif iSubplot == 13
        % Center subplot 13 in the first half of the last row
            ax(iSubplot).pos(1) = 0.20 - ax(iSubplot).pos(3)/2;
        elseif iSubplot == 14
        % Center subplot 14 in the second half of the last row
            ax(iSubplot).pos(1) = 0.75 - ax(iSubplot).pos(3)/2;
        end
    
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColourMap,...
            caxisMin(iSubplot),caxisMax(iSubplot));
        title(labelVars(iSubplot),'FontSize',12) 
        
        cb = colorbar('Location','eastoutside');   
        cb.Position(1) = cb.Position(1) + 0.055;
        cb.Position(3) = 0.015; % width
        cb.Position(4) = 0.099; % length
        
        cbStringModified = strrep(cbString(iSubplot), '(log10)', ''); % remove 'log10' from the string
        cb.Label.String = cbStringModified; % set colorbar label without 'log10'
        
        % Center the label above the vertical colorbar
        cb.Label.Rotation = 0;
        cb.Label.HorizontalAlignment = 'center';    
        cb.Label.VerticalAlignment = 'bottom'; 
        cb.Label.Position(1) = 0.60;
        cb.Label.Position(2) = caxisMax(iSubplot); 

    end % iSubplot

    saveFigure(figureName)

elseif (nSubplots == 15 && ~isCommonColourBar)

    figure()
    set(gcf,'Units','Normalized','Position',[0.0 0.0 0.50 0.89],'Color','w') 
    haxis = zeros(nSubplots,1);

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(5,3,iSubplot,'Spacing',0.015,'Padding',0.015,'Margin',0.027);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1 || iSubplot == 4 || iSubplot == 7 || iSubplot == 10 || iSubplot == 13)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.03; 
        elseif (iSubplot == 2 || iSubplot == 5 || iSubplot == 8 || iSubplot == 11 || iSubplot == 14)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02; 
        elseif (iSubplot == 3 || iSubplot == 6 || iSubplot == 9 || iSubplot == 12 || iSubplot == 15)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.01;
        end
    
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

        mydata = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lonVector,latVector,mydata,myColourMap,...
            caxisMin(iSubplot),caxisMax(iSubplot));
        title(labelVars(iSubplot),'FontSize',12) 
        
        cb = colorbar('Location','eastoutside');   
        cb.Position(1) = cb.Position(1) + 0.055;
        cb.Position(3) = 0.015; % width
        cb.Position(4) = 0.099; % length
        
        cbStringModified = strrep(cbString(iSubplot), '(log10)', ''); % remove 'log10' from the string
        cb.Label.String = cbStringModified; % set colorbar label without 'log10'
        
        % Center the label above the vertical colorbar
        cb.Label.Rotation = 0;
        cb.Label.HorizontalAlignment = 'center';    
        cb.Label.VerticalAlignment = 'bottom'; 
        cb.Label.Position(1) = 0.60;
        cb.Label.Position(2) = caxisMax(iSubplot); 

    end % iSubplot

    saveFigure(figureName)
    
elseif (nSubplots == 18)

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.55 0.85],'Color','w')
    haxis = zeros(nSubplots,1);

    [oceanVar,caxisMin,caxisMax,ticks,tickLabels] = assessLogarithmicTransformation(oceanVar);
    
    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(6,3,iSubplot,'Spacing',0.005,'Padding',0.015,'Margin',0.04);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        ax(iSubplot).pos(2) =  ax(iSubplot).pos(2) + 0.03;
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

        myData = oceanVar(:,:,iSubplot);
        plotPcolorMap(haxis(iSubplot),lon,lat,myData,myColourmap,caxisMin,caxisMax);
        title(labelVars(iSubplot),'FontSize',11)     

    end

    cb = colorbar(haxis(nSubplots));
    cb.Location = 'southoutside';
    cb.Position(3) = 0.40; % length
    cb.Position(4) = 0.01; % width
    cb.Position(1) = 0.50-cb.Position(3)/2;
    cb.Position(2) = 0.05; 
    cb.Label.String = cbString;
    
    % Set axis limits and ticks
    caxis([caxisMin caxisMax]);
    cb.Ticks = ticks;
    cb.TickLabels = tickLabels;
    cb.FontSize = 8; % font size of tick labels
    cb.Label.FontSize = 10; % font size of colour bar string

    saveFigure(figureName)

end

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function [outputData,caxisMin,caxisMax,ticks,tickLabels] =...
    assessLogarithmicTransformation(originalData)

% originalData = oceanVar;

    % This function:
    % - handling of logarithmic and linear scaling separately
    % - dynamically adjusts step sizes to ensure a maximum of 8 ticks
    % - handles small ranges by refining limits and step sizes
    % - robust percentile handling to prevent extreme outliers from skewing the limits

    useLogScale = false;
    
    originalData(originalData==0) = NaN; % make sure there are no zeros(log10 cannot cope with those)

    % Compute log range and skewness
    logRange = log10(max(originalData(:))) - log10(min(originalData(:)));
    dataSkewness = skewness(originalData(:));
    % figure()    
    % histogram(originalData, 20);     % Plot histogram with 20 bins
    % xlabel('Data values');
    % ylabel('Frequency');
    % title('Histogram of Data')

    % Use log scale if data spans at least 4 orders of magnitude OR is highly skewed
    if logRange > 3 || dataSkewness > 1
        useLogScale = true;
    end

    % **LOGARITHMIC SCALE HANDLING**
    if useLogScale
        outputData = log10(originalData);
        logData = log10(originalData(:));

        % Compute 5th and 95th percentiles to exclude extreme outliers
        logMin = prctile(logData, 5);
        logMax = prctile(logData, 95);
          
        % Round to nearest power of 10
        caxisMin = floor(logMin);
        caxisMax = floor(logMax);

        % Adjust max limit slightly if needed
        if (max(logData) - caxisMax) < 0.5
            caxisMax = caxisMax + 1;
        end

        % Determine step size to ensure max 8 ticks
        nTicks = caxisMax - caxisMin + 1;
        if nTicks > 8
            stepSize = ceil((caxisMax - caxisMin) / 8);
        else
            stepSize = 1;
        end

        % Generate tick positions
        ticks = caxisMin:stepSize:caxisMax;

        % Create labels
        tickLabels = arrayfun(@(x) sprintf('10^{%d}',x), ticks, 'UniformOutput', false);
        
    else
        % **LINEAR SCALE HANDLING**
        outputData = originalData;

        % Compute 5th and 95th percentiles for robust min/max
        linearMin = prctile(originalData(:), 5);
        linearMax = prctile(originalData(:), 95);
        
        % Ensure min and max are not identical
        if linearMax == linearMin
            linearMax = linearMax * 1.1;
        end
        
        % **Handle very small or very large numbers**
        range = linearMax - linearMin;
        if range < 1e-6
            
            % Define caxisMin and caxisMax properly in linear space
            caxisMin = floor(linearMin * 10^ceil(-log10(linearMax))) / 10^ceil(-log10(linearMax));
            caxisMax = ceil(linearMax * 10^ceil(-log10(linearMax))) / 10^ceil(-log10(linearMax));

            % Ensure at least 3 ticks
            nTicks = max(3, min(8, ceil(log10(caxisMax) - log10(caxisMin)) + 1));

            % Generate tick positions with appropriate spacing
            ticks = linspace(caxisMin, caxisMax, nTicks);

            % Format tick labels using scientific notation
            tickLabels = arrayfun(@(x) sprintf('%.1e', x), ticks, 'UniformOutput', false);
    
        else
    
            % Use standard rounding for normal ranges
            roundFactor = 10^floor(log10(range));
            caxisMin = floor(linearMin / roundFactor) * roundFactor;
            caxisMax = ceil(linearMax / roundFactor) * roundFactor;

            % Ensure max > min
            if caxisMax <= caxisMin
                caxisMax = caxisMin + roundFactor;
            end

            % Determine step size dynamically (max 8 ticks)
            stepSizes = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
            stepSize = stepSizes(find((caxisMax - caxisMin) ./ stepSizes <= 8, 1, 'first'));

            % Generate tick positions with at least 3 ticks but no more
            % than 8
            ticks = caxisMin:stepSize:caxisMax;
            if numel(ticks) < 3
                ticks = linspace(caxisMin, caxisMax, 3);
            elseif numel(ticks) > 8
                ticks = linspace(caxisMin, caxisMax, 8);
            end
            
            % Format tick labels
            % Ensure caxisMin is not zero (log10(0) is undefined)
            if caxisMin == 0
                caxisMin = min(ticks(ticks > 0)); % Set to the smallest positive tick
            end
            if abs(log10(caxisMin)) >= 4 || abs(log10(caxisMax)) >= 4
                tickLabels = arrayfun(@(x) sprintf('%.1e', x), ticks, 'UniformOutput', false);
            else
                tickLabels = arrayfun(@(x) sprintf('%g', x), ticks, 'UniformOutput', false);
            end

        end

    end

end % assessLogarithmicTransformation

end % plotOceanVariableMaps
