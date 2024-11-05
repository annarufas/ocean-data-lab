
% ======================================================================= %
%                                                                         %
% This script computes and visually compares the annual mean of various   % 
% ocean variables from multiple data sources, facilitating cross-source   %
% comparison of spatial and temporal variations.                          %
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

fullpathInputDataDir = fullfile('data','processed');

% Chla
fullpathChlaAquaModisFile = fullfile(fullpathInputDataDir,'chla_aquamodis.mat');
fullpathChlaCmemsFile     = fullfile(fullpathInputDataDir,'chla_cmems_bgc.mat');
fullpathChlaOccciFile     = fullfile(fullpathInputDataDir,'chla_occci.mat');

% NPP
fullpathNppBicepFile       = fullfile(fullpathInputDataDir,'npp_bicep.mat');
fullpathNppCmemsFile       = fullfile(fullpathInputDataDir,'npp_cmems_bgc.mat');
fullpathNppCafeModisFile   = fullfile(fullpathInputDataDir,'npp_cafe_modis.mat');
fullpathNppCafeSeawifsFile = fullfile(fullpathInputDataDir,'npp_cafe_seawifs.mat');
fullpathNppCbpmFile        = fullfile(fullpathInputDataDir,'npp_cbpm_modis.mat');
fullpathNppVgpmFile        = fullfile(fullpathInputDataDir,'npp_vgpm_modis.mat');

% MLD
fullpathMldIfremerFile = fullfile(fullpathInputDataDir,'mld_ifremer.mat');
fullpathMldCmemsFile   = fullfile(fullpathInputDataDir,'mld_cmems_phys.mat');

% Seawater temperature
fullpathTempCmemsFile = fullfile(fullpathInputDataDir,'temp_cmems_phys.mat');
fullpathTempWoaFile   = fullfile(fullpathInputDataDir,'temp_monthly_woa23.mat');

% PAR0
fullpathPar0AquaModisFile  = fullfile(fullpathInputDataDir,'par0_aquamodis.mat');
fullpathPar0CalculatedFile = fullfile(fullpathInputDataDir,'par0_monthly_calculated.mat');

% Colour maps
myColorMapVar = flipud(brewermap(100, 'RdYlBu')); % flip to have blue colours for low values
myColorMapDiff = flipud(brewermap(100, 'RdBu')); % flip to have blue colours for low values

% Grid used for regridding data to common resolution
fullpathGridFile = fullfile('data','raw','grid','grid_GEBCO_2160_1080.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD THE GRID AND PREPARE INTERPOLANTS
% -------------------------------------------------------------------------

load(fullpathGridFile,'Xbb','Ybb','x','y','ixBb','iyBb') 
qLons = double(x);
qLats = double(y);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CHLA COMPARISON
% -------------------------------------------------------------------------

filePaths = {fullpathChlaAquaModisFile,fullpathChlaCmemsFile,fullpathChlaOccciFile};
suffixes = {'cmems','occci','aquamodis'}; % matches those in 'filePaths'

% Extract data from the various sources listed in 'filePaths'
varStack = [];
titleStrStack = cell(1, length(filePaths));
for i = 1:length(filePaths)
    % Load the data
    data = load(filePaths{i},'chla','chla_lon','chla_lat'); % adjust accordingly
    % Regrid
    Dout = regridToCommonResolution(mean(data.chla,3,'omitnan'),...
        data.chla_lon,data.chla_lat,qLons,qLats);
    % Outputs
    titleStrStack{i} = upper(suffixes{i});
    varStack(:,:,i) = Dout;
    assignin('base', ['chla_' suffixes{i}], Dout);
end

% Compute chla differences
chla_comp1 = chla_cmems - chla_occci;
chla_comp2 = chla_cmems - chla_aquamodis;
chla_comp3 = chla_occci - chla_aquamodis;
allVars = cat(3,varStack,chla_comp1,chla_comp2,chla_comp3);

% Plot subplots
titleStr = [titleStrStack,'[CMEMS] - [OC-CCI]','[CMEMS] - [Aqua-MODIS]','[OC-CCI] - [Aqua-MODIS]'];
cbString = {'Annual mean chla (mg m^{-3})','Difference'};
organiseSubplotArrangement(size(allVars,3),allVars,qLons,qLats,...
    myColorMapVar,myColorMapDiff,[0, -0.1],[1, 0.1],titleStr,cbString)
saveFigure([],'fig_comparison_chla')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - NPP COMPARISON
% -------------------------------------------------------------------------

filePaths = {fullpathNppCmemsFile,fullpathNppBicepFile,fullpathNppCafeModisFile};
suffixes = {'cmems','bicep','cafe'}; % matches those in 'filePaths'

% Extract data from the various sources listed in 'filePaths'
varStack = [];
titleStrStack = cell(1, length(filePaths));
for i = 1:length(filePaths)    
    % Load the data
    data = load(filePaths{i},'npp_avg','npp_lon','npp_lat'); % adjust accordingly
    % Regrid
    Dout = regridToCommonResolution(mean(data.npp_avg,3,'omitnan'),...
        data.npp_lon,data.npp_lat,qLons,qLats);   
    % Outputs
    titleStrStack{i} = upper(suffixes{i});
    varStack(:,:,i) = Dout;
    assignin('base', ['npp_' suffixes{i}], Dout);
end

% Compute differences
npp_comp1 = npp_cmems - npp_bicep;
npp_comp2 = npp_cmems - npp_cafe; 
npp_comp3 = npp_bicep - npp_cafe;
allVars = cat(3,varStack,npp_comp1,npp_comp2,npp_comp3);

% Plot subplots
titleStr = [titleStrStack,'[CMEMS] - [BICEP]','[CMEMS] - [CAFE]','[BICEP] - [CAFE]'];
cbString = {'Annual mean NPP (mg C m^{-2})','Difference'};
organiseSubplotArrangement(size(allVars,3),allVars,qLons,qLats,...
    myColorMapVar,myColorMapDiff,[0, -300],[1000, 300],titleStr,cbString)
saveFigure([],'fig_comparison_npp')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - MLD COMPARISON
% -------------------------------------------------------------------------

filePaths = {fullpathMldIfremerFile,fullpathMldCmemsFile};
suffixes = {'cmems','ifremer'}; % matches those in 'filePaths'

% Extract data from the various sources listed in 'filePaths'
varStack = [];
titleStrStack = cell(1, length(filePaths));
for i = 1:length(filePaths)
    % Load the data
    data = load(filePaths{i},'mld','mld_lon','mld_lat'); % adjust accordingly
    % Regrid
    Dout = regridToCommonResolution(mean(data.mld,3,'omitnan'),...
        data.mld_lon,data.mld_lat,qLons,qLats);
    % Outputs
    titleStrStack{i} = upper(suffixes{i});
    varStack(:,:,i) = Dout;
    assignin('base', ['mld_' suffixes{i}], Dout);
end

% Compute differences
mld_comp = mld_cmems - mld_ifremer;
allVars = cat(3,varStack,mld_comp);

% Plot subplots
titleStr = [titleStrStack,'[CMEMS] - [IFREMER]'];
cbString = {'Annual mean MLD (m)','Difference'};
organiseSubplotArrangement(size(allVars,3),allVars,qLons,qLats,...
    myColorMapVar,myColorMapDiff,[20, -30],[200, 30],titleStr,cbString)
saveFigure([],'fig_comparison_mld')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - SEAWATER TEMPERATURE COMPARISON
% -------------------------------------------------------------------------

filePaths = {fullpathTempCmemsFile,fullpathTempWoaFile};
suffixes = {'cmems','woa'}; % matches those in 'filePaths'
targetDepth = 50; % m

% Extract data from the various sources listed in 'filePaths'
varStack = [];
titleStrStack = cell(1, length(filePaths));
for i = 1:length(filePaths)
    % Load the data (adjust accordingly)
    if (i == 1)
        data = load(filePaths{i},'temp','temp_lon','temp_lat','temp_depth');
        [~, idxDepth] = min(abs(data.temp_depth - targetDepth));
        lon = data.temp_lon;
        lat = data.temp_lat;
    elseif (i == 2)
        data = load(filePaths{i},'temp','woa_lon','woa_lat','woa_depth_temp');
        [~, idxDepth] = min(abs(data.woa_depth_temp - targetDepth));
        lon = data.woa_lon;
        lat = data.woa_lat;
    end
    % Regrid
    selectSlice = squeeze(data.temp(:,:,idxDepth,:));
    Dout = regridToCommonResolution(mean(selectSlice,3,'omitnan'),...
        lon,lat,qLons,qLats);
    % Outputs
    titleStrStack{i} = upper(suffixes{i});
    varStack(:,:,i) = Dout;
    assignin('base', ['temp_' suffixes{i}], Dout);
end

% Compute differences
temp_comp = temp_cmems - temp_woa;
allVars = cat(3,varStack,temp_comp);

% Plot subplots
titleStr = [titleStrStack,'[CMEMS] - [WOA]'];
cbString = {'Annual mean temperature at 50 m (m)','Difference'};
organiseSubplotArrangement(size(allVars,3),allVars,qLons,qLats,...
    myColorMapVar,myColorMapDiff,[0, -2],[25, 2],titleStr,cbString)
saveFigure([],'fig_comparison_temp')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - PAR0 COMPARISON
% -------------------------------------------------------------------------

filePaths = {fullpathPar0AquaModisFile,fullpathPar0CalculatedFile};
suffixes = {'aquamodis','calculated'}; % matches those in 'filePaths'

% Extract data from the various sources listed in 'filePaths'
varStack = [];
titleStrStack = cell(1, length(filePaths));
for i = 1:length(filePaths)
    % Load the data (adjust accordingly)
    if (i == 1)
        data = load(filePaths{i},'par0','par0_lon','par0_lat');
        dataVar = data.par0;
    elseif (i == 2)
        data = load(filePaths{i},'par0clim','par0_lon','par0_lat');
        dataVar = data.par0clim;
    end
    % Regrid
    Dout = regridToCommonResolution(mean(dataVar,3,'omitnan'),...
        data.par0_lon,data.par0_lat,qLons,qLats);
    % Outputs
    titleStrStack{i} = upper(suffixes{i});
    varStack(:,:,i) = Dout;
    assignin('base', ['par0_' suffixes{i}], Dout);
end

% Compute differences
par0_comp = par0_aquamodis - par0_calculated;
allVars = cat(3,varStack,par0_comp);

% Plot subplots
titleStr = [titleStrStack,'[Aqua-MODIS] - [calculated]'];
cbString = {'Annual mean PAR0 (W m^{-2})','Difference'};
organiseSubplotArrangement(size(allVars,3),allVars,qLons,qLats,...
    myColorMapVar,myColorMapDiff,[50, -30],[150, 30],titleStr,cbString)
saveFigure([],'fig_comparison_par0')

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function [Dout] = regridToCommonResolution(annualMean,lon,lat,qLons,qLats)
    
    % Determine dimension order
    [~,iDimLat] = min(size(annualMean));
    [~,iDimLon] = max(size(annualMean));
    isLonFirstDim = (iDimLon == 1);
    
    % Regrid to a common resolution 
    if isLonFirstDim
        [X,Y] = ndgrid(lon,lat); % original array
        [qX,qY] = ndgrid(qLons,qLats); % query points for interpolation 
    else
        [X,Y] = ndgrid(lat,lon); % original array
        [qX,qY] = ndgrid(qLats,qLons); % query points for interpolation 
    end
    F = griddedInterpolant(X,Y,annualMean);
    qAnnualMean = F(qX,qY); 
    
    % Permute dimensions so all output arrays are lat x lon
    Dperm = permute(qAnnualMean,[iDimLat iDimLon]);
    Dout = Dperm;
    
end % regridToCommonResolution

% *************************************************************************

function organiseSubplotArrangement(nSubplots,oceanVar,qLons,qLats,...
    myColorMapVar,myColorMapDiff,caxisMin,caxisMax,titleStr,cbString)

    if (nSubplots == 6)

        figure()
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.50 0.48],'Color','w')
        haxis = zeros(nSubplots,1);

        for iSubplot = 1:nSubplots

            haxis(iSubplot) = subaxis(2,3,iSubplot,'Spacing',0.00,'Padding',0.005,'Margin',0.00);
            ax(iSubplot).pos = get(haxis(iSubplot),'Position');
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.02;
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

            if (iSubplot < 4)
                plotPcolorMap(haxis(iSubplot),qLons,qLats,oceanVar(:,:,iSubplot),...
                    myColorMapVar,caxisMin(1),caxisMax(1));
                title(titleStr(iSubplot),'FontSize',13) 
                if (iSubplot == 2)
                    cb = colorbar(haxis(2));
                    cb.Location = 'southoutside';
                    cb.Position(2) = cb.Position(2) - 0.02;  
                    cb.Label.String = cbString(1);
                    cb.FontSize = 12;
                end
            elseif (iSubplot >= 4)
                plotPcolorMap(haxis(iSubplot),qLons,qLats,oceanVar(:,:,iSubplot),...
                    myColorMapDiff,caxisMin(2),caxisMax(2));
                title(titleStr(iSubplot),'FontSize',13) 
                if (iSubplot == 5)
                    cb = colorbar(haxis(5));
                    cb.Location = 'southoutside';
                    cb.Position(2) = cb.Position(2) - 0.02;  
                    cb.Label.String = cbString(2);
                    cb.FontSize = 12;
                end
            end

        end

    elseif (nSubplots == 3)

        figure()
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.48],'Color','w')
        haxis = zeros(nSubplots,1);

        for iSubplot = 1:nSubplots

            haxis(iSubplot) = subaxis(2,2,iSubplot,'Spacing',0.00,'Padding',0.005,'Margin',0.00);
            ax(iSubplot).pos = get(haxis(iSubplot),'Position');
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.04;
            if (iSubplot == 3) % move 3rd subplot to the middle of the 2nd row
                ax(iSubplot).pos(1) = 0.5 - 0.25; % center in the 2nd row
            end
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

            if (iSubplot < 3)
                plotPcolorMap(haxis(iSubplot),qLons,qLats,oceanVar(:,:,iSubplot),...
                    myColorMapVar,caxisMin(1),caxisMax(1));
                title(titleStr(iSubplot),'FontSize',13) 
                if (iSubplot == 1)
                    cb = colorbar(haxis(1));
                    cb.Location = 'southoutside';
                    cb.Position(1) = 0.25; % center it between subplot 1 and 2
                    cb.Position(2) = cb.Position(2) - 0.02;  
                    cb.Label.String = cbString(1);
                    cb.FontSize = 12;
                end
            elseif (iSubplot >= 3)
                plotPcolorMap(haxis(iSubplot),qLons,qLats,oceanVar(:,:,iSubplot),...
                    myColorMapDiff,caxisMin(2),caxisMax(2));
                title(titleStr(iSubplot),'FontSize',13) 
                if (iSubplot == 3)
                    cb = colorbar(haxis(3));
                    cb.Location = 'southoutside';
                    cb.Position(1) = 0.25; % center it under subplot 3
                    cb.Position(2) = cb.Position(2) - 0.02;  
                    cb.Label.String = cbString(2);
                    cb.FontSize = 12;
                end
            end

        end

    end

end % organiseSubplotArrangement

% *************************************************************************
