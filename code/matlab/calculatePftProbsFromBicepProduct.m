
% ======================================================================= %
%                                                                         %
% This script processes BICEP Cphyto data (mg C m-3), monthly global      %
% arrays for 1998-2020 and calculates PFT seeding probabilities that can  %
% be used to force SLAMS. BICEP uses this size classification (after      %
% Brewin et al. (2015):                                                   %
%   Microphytoplankton: > 20 um in diameter                               %
%   Nanophytoplankton: 2-20 um in diameter                                %
%   Picophytoplankton: < 2 um in diameter                                 %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 14 August 2024                                %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./modelresources/external/')); 
addpath(genpath('./code/matlab/')); 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Path to dataset
fullpathBicepCphytoDataDir = '/Users/Anna/LocalDocuments/Academic/Projects/ocean-data-lab/data/raw/BICEP/BICEP_Cphyto/';
figuresDir = 'bicep';

% Time vector
yearsVector = 1998:2020; 

% Variabels to extract
oceanColourVars = {'C_phyto','C_microphyto','C_nanophyto','C_picophyto'};

% Plotting properties
mycolormap = [ones(1,3); jet(1000)];

labelFourPfts = {'{Diatoms}',...
                 '{Flagellated phytoplankton}',...
                 '{Coccolithophores}',...
                 '{Picophytoplankton}'};

labelThreePfts = {'{Microphytoplankton}',...
                  '{Nanophytoplankton}',...
                  '{Picophytoplankton}'};
              
labelMonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT BICEP DATA
% -------------------------------------------------------------------------

for iMonth = 1:12
    disp(iMonth)
    for iYear = 1:numel(yearsVector)
        
        thisYear = yearsVector(iYear);
        yearFolderPath = fullfile(fullpathBicepCphytoDataDir, num2str(thisYear));
        fileNames = dir(fullfile(yearFolderPath,'*.nc'));

        filePath = fullfile(yearFolderPath, fileNames(iMonth).name);
        %S = ncinfo(filePath); % short summary

        % Read longitude and latitude
        % Note: latitude might be unsorted
        lat = ncread(filePath,'latitude');
        lon = ncread(filePath,'longitude');
        
        if (iYear == 1 && iMonth == 1)
            pftClimatology = NaN(numel(lat),numel(lon),12,numel(oceanColourVars));    
        end

        % Read data and permute lat and lon
        D = zeros(numel(lat),numel(lon),numel(oceanColourVars));
        for iVar = 1:numel(oceanColourVars)
            Dtmp = ncread(filePath,oceanColourVars{iVar});
            Dperm = permute(Dtmp,[2 1]); % swap lat and lon
            D(:,:,iVar) = Dperm;
        end

        % Sort lat and lon to have monotonically increasing values and 
        % apply the sorting to the data
        [lat_sort, sortLatIdx] = sort(lat);
        [lon_sort, sortLonIdx] = sort(lon);
        D_sort = D(sortLatIdx,sortLonIdx,:); % lat x lon x vars

        % Initialise the output data array on the first iteration
        if (iYear == 1)    
            D_out = NaN(numel(lat_sort),numel(lon_sort),numel(yearsVector),numel(oceanColourVars),'single');
        end

        % Store the sorted data in the output array at the calculated time index
        D_out(:,:,iYear,:) = D_sort(:,:,:);

    end % iYear
    
    pftClimatology(:,:,iMonth,:) = mean(D_out(:,:,:,:),3,'omitnan'); % mg C m-3
   
end % iMonth

% Save the data
pft_bicep_lon = lon;
pft_bicep_lat = lat;
pft_bicep = pftClimatology;
save(fullfile('.','data','raw','pft_bicep.mat'),...
    'pft_bicep','pft_bicep_lon','pft_bicep_lat','-v7.3')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE BICEP RELATIVE ABUNDANCES
% -------------------------------------------------------------------------

load(fullfile('.','data','raw','pft_bicep.mat'),'pft_bicep','pft_bicep_lon','pft_bicep_lat')

pftBiomassMonthly = pft_bicep(:,:,:,2:4); % the 1st group is total phytoplankton biomass; 2, 3, 4 are Cmicro, Cnano and Cpico
pftBiomassAnnual = squeeze(mean(pftBiomassMonthly,3,'omitnan')); % mg C m-3
pftBiomassAnnualTotal = sum(pftBiomassAnnual,3,'omitnan'); % mg C m-3

% Calculate relative biomass where total biomass is greater than 0
pftBiomassRelative = zeros(numel(pft_bicep_lat),numel(pft_bicep_lon),3);
for iPart = 1:3
    pftBiomassRelative(:,:,iPart) = 100.*pftBiomassAnnual(:,:,iPart)./pftBiomassAnnualTotal(:,:);
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - PLOT BICEP DATA
% -------------------------------------------------------------------------

longrid = pft_bicep_lon;
latgrid = pft_bicep_lat;

for iFigure = 1:2
    
    if (iFigure == 1)
        
        oceanVar = pftBiomassAnnual;
        caxismin = 0;
        caxismax = 15;
        cbString = 'Phytoplankton standing stock over z_{eu} (mg C m^{-3})';
        figureName = 'bicep_three_phyto_biomass';
        
        plotOceanVariableMaps(oceanVar,longrid,latgrid,mycolormap,cbString,...
            caxismin,caxismax,labelThreePfts,figuresDir,figureName,[])
    
    elseif (iFigure == 2)

        oceanVar = pftBiomassRelative;
        caxismin = 0;
        caxismax = 50;
        cbString = 'Relative biomass (%)';
        figureName = 'bicep_three_phyto_biomass_relative';
        
        plotOceanVariableMaps(oceanVar,longrid,latgrid,mycolormap,cbString,...
            caxismin,caxismax,labelThreePfts,figuresDir,figureName,[])

    end

end % iFigure


% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - IN BICEP, REPARTITION MICROPHYTO INTO DIATOMS AND FLAGELLATES
% -------------------------------------------------------------------------

% For BICEP, let's assume a 1/3 of the microphytoplankton carbon is large
% flagellates and the rest is diatoms

factorDiatoms = 2/3;

pftBiomassMonthlyDiatoms = pftBiomassMonthly(:,:,:,1).*factorDiatoms;
pftBiomassMonthlyLargeFlagel = pftBiomassMonthly(:,:,:,1).*(1-factorDiatoms);

% Assign the calculated and original PFTs to the new array
fourPftBiomassMonthly = zeros(numel(pft_bicep_lat),numel(pft_bicep_lon),12,4);
fourPftBiomassMonthly(:,:,:,1) = pftBiomassMonthlyDiatoms;
fourPftBiomassMonthly(:,:,:,2) = pftBiomassMonthlyLargeFlagel;
fourPftBiomassMonthly(:,:,:,3) = pftBiomassMonthly(:,:,:,2);
fourPftBiomassMonthly(:,:,:,4) = pftBiomassMonthly(:,:,:,3);

% Do calculations
fourPftBiomassAnnual = squeeze(mean(fourPftBiomassMonthly,3,'omitnan')); % mg C m-3
fourPftBiomassAnnualTotal = sum(fourPftBiomassAnnual,3,'omitnan');
fourPftBiomassRelative = zeros(numel(pft_bicep_lat),numel(pft_bicep_lon),4);
for iPart = 1:4
    fourPftBiomassRelative(:,:,iPart) = 100.*fourPftBiomassAnnual(:,:,iPart)./fourPftBiomassAnnualTotal(:,:);
end

% Plot

mycolorbar = [ones(1,3); jet(1000)];

labelSubplots = {'{Diatoms}',...
                 '{Flagellated phytoplankton}',...
                 '{Coccolithophores}',...
                 '{Picophytoplankton}'};
               
nSubplots = length(labelSubplots);

for iFigure = 1:2
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.45],'Color','w')
    haxis = zeros(nSubplots,1);
    
    if (iFigure == 1)
        mydata = fourPftBiomassAnnual;
        caxismin = 0;
        caxismax = 15;
        cbString = 'Biomass in z_{eu} (mg C m^{-3})';
        figName = 'bicep_four_phyto_biomass';
    elseif (iFigure == 2)
        mydata = fourPftBiomassRelative;
        caxismin = 0;
        caxismax = 50;
        cbString = 'Relative biomass (%)';
        figName = 'bicep_four_phyto_biomass_relative';
    end

    for iSubplot = 1:nSubplots

        haxis(iSubplot) = subaxis(2,2,iSubplot,'Spacing',0.01,'Padding',0.02,'Margin',0.04);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.05;
        elseif (iSubplot == 2)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02; 
        elseif (iSubplot == 3)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.05;
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.08;
        elseif (iSubplot == 4)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02; 
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.08;
        end
        set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 

        m_proj('Robinson','clongitude',-63);

        m_pcolor(pft_bicep_lon,pft_bicep_lat,mydata(:,:,iSubplot))
        if (sum(sum(pft_bicep_lon < 0)) > 0)
            hold on
            m_pcolor(pft_bicep_lon-360,pft_bicep_lat,mydata(:,:,iSubplot))
        end
        caxis([caxismin caxismax])

        shading flat; colormap(mycolorbar)
        m_grid('xticklabels',[],'yticklabels',[],'ticklen',0,'linestyle','none'); 
        m_coast('patch',[.7 .7 .7],'edgecolor','black');

        title(labelSubplots(iSubplot),'FontSize',13)

    end
    
    cb = colorbar(haxis(4));
    cb.Location = 'southoutside';
    cb.Position(1) = 0.50-cb.Position(3)/2; 
    cb.Position(2) = ax(3).pos(2) - cb.Position(4);
    cb.Label.String = cbString;
    cb.FontSize = 12;

    exportgraphics(gcf,fullfile(figuresDir,strcat(figName,'.png')),'Resolution',600)

end % iFigure