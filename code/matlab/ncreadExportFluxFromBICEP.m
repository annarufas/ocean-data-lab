
% ======================================================================= %
%                                                                         %
%            POC export flux climatology from BICEP project               %
%                                                                         %
% This script reads in POC export flux data from the Biological Pump and  %
% Carbon Exchange Processes (BICEP) project.                              %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://data.ceda.ac.uk/neodc/bicep/data/oceanic_export_production/ %
%   - Time resolution: monthly climatology                                %
%   - Space resolution: 9 km resolution (2160 x 4320)                     %
%   - Version: v4.2, 1998–2019                                            %
%   - Units: mg C m-2 d-1                                                 %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 28 May 2025                                   %
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

fullpathInputExFluxDir = fullfile('data','raw','BICEP','BICEP_export_flux');
fullpathOutputExFluxDunneFile = fullfile('data','processed','ef_dunne_bicep.mat');
fullpathOutputExFluxHensonFile = fullfile('data','processed','ef_henson_bicep.mat');
fullpathOutputExFluxLiFile = fullfile('data','processed','ef_li_bicep.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

% Time vector
yearsVector = 1998:2019;

for iMonth = 1:12
    for iYear = 1:numel(yearsVector)       
        thisYear = yearsVector(iYear);
        yearFolderPath = fullfile(fullpathInputExFluxDir, num2str(thisYear));
        fileNames = dir(fullfile(yearFolderPath,'*.nc'));

        filePath = fullfile(yearFolderPath, fileNames(iMonth).name);
        %S = ncinfo(filePath); % short summary

        % Read longitude and latitude
        lat = ncread(filePath,'lat');
        lon = ncread(filePath,'lon');
        
        if (iYear == 1 && iMonth == 1)
            exportflux = NaN(numel(lat),numel(lon),3,12);
        end

        % Read data and permute lat and lon
        modelTags = {'EP_Dunne','EP_Henson','EP_Li'};
        D = zeros(numel(lat),numel(lon),numel(modelTags));
        for iModel = 1:numel(modelTags)
            Dtmp = ncread(filePath, modelTags{iModel});
            D(:,:,iModel) = permute(Dtmp, [2 1]); % swap lat and lon
        end

        % Sort lat and lon to have monotonically increasing values and 
        % apply the sorting to the data
        [lat_sort, sortLatIdx] = sort(lat);
        [lon_sort, sortLonIdx] = sort(lon);
        D_sort = D(sortLatIdx,sortLonIdx,:); % lat x lon x vars

        % Initialise the output data array on the first iteration
        if (iYear == 1)    
            D_out = NaN(numel(lat_sort),numel(lon_sort),numel(modelTags),numel(yearsVector),'single');
        end

        % Store the sorted data in the output array at the calculated time index
        D_out(:,:,:,iYear) = D_sort(:,:,:);

    end % iYear
    
    exportflux(:,:,:,iMonth) = mean(D_out(:,:,:,:),4,'omitnan'); % mg C m-2 d-1
   
end % iMonth

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CHECK VALUES AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(exportflux(:),100);

% Save the data
ef_lon = lon_sort;
ef_lat = lat_sort;
ef_dunne_avg = squeeze(exportflux(:,:,1,:));
ef_henson_avg = squeeze(exportflux(:,:,2,:));
ef_li_avg = squeeze(exportflux(:,:,3,:));
save(fullpathOutputExFluxDunneFile,'ef_dunne_avg','ef_lat','ef_lon','-v7.3')
save(fullpathOutputExFluxHensonFile,'ef_henson_avg','ef_lat','ef_lon','-v7.3')
save(fullpathOutputExFluxLiFile,'ef_li_avg','ef_lat','ef_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputExFluxDunneFile,[],'mg C m^{-2} d^{-1}',...
    0,150,true,'fig_monthly_ef_dunne_bicep','POC export flux, Dunne, BICEP')
prepareDataForPlotting(fullpathOutputExFluxHensonFile,[],'mg C m^{-2} d^{-1}',...
    0,150,true,'fig_monthly_ef_henson_bicep','POC export flux, Henson, BICEP')
prepareDataForPlotting(fullpathOutputExFluxLiFile,[],'mg C m^{-2} d^{-1}',...
    0,150,true,'fig_monthly_ef_li_bicep','POC export flux, Li, BICEP')

%% Plot annual

    myColourMap = jet(1000);
    caxisMin = 0; % mg C m-2 d-2
    caxisMax = 150; % mg C m-2 d-2
    cbString = 'POC export flux, Li, BICEP (mg C m^{-2} d^{-1})';


    % The monthly NPP data arrays will be regridded to a commmon 1º lat x 1º lon
    qLat = linspace(-90,90,180); % query points for interpolation
    qLon = linspace(-180,180,360); % query points for interpolation
    [qX, qY, qT] = ndgrid(qLat, qLon, 1:12); % query grid with 12 time steps

    % Output annual arrays
    nppModelledAnnualMeanAll = NaN(numel(qLat),numel(qLon),nNppModels*nGapFillingMethods);
    nppModelledAnnualMeanGapFilled = NaN(numel(qLat),numel(qLon),nNppModels);
    globalNppStockSummaryGapFilled = cell(nNppModels,1);

    counter = 1;
    for iModel = 1:nNppModels 
        for iMethod = 1:nGapFillingMethods
            switch iMethod
                case 1, label = 'plain';
                case 2, label = 'gapfilled';
            end

            % Get data from the structure
            fileName = erase(filenameModelledNpp{iModel},'.mat');
            fieldName = [fileName, '_', label]; 
            data = nppModelClimatologyStruct.(fieldName).data;
            lat = nppModelClimatologyStruct.(fieldName).lat;
            lon = nppModelClimatologyStruct.(fieldName).lon;

            % Regrid to common 1º lat x 1º lon
            [Xpp, Ypp, Tpp] = ndgrid(lat, lon, (1:12)'); % original grid for the current dataset
            Favg = griddedInterpolant(Xpp, Ypp, Tpp, data);
            qNppAvg = squeeze(Favg(qX, qY, qT));

            % Compute the annual mean and store in the output array
            nppModelledAnnualMeanAll(:,:,counter) = mean(qNppAvg,3,'omitnan');
            counter = counter + 1;
            
            if (iMethod == 2)
                iDataset = (iModel - 1)*nGapFillingMethods + iMethod;
                titleStr = strrep(globalNppStockSummary{iDataset}, ' gapfilled', '');
                nppModelledAnnualMeanGapFilled(:,:,iModel) = mean(qNppAvg,3,'omitnan');
                globalNppStockSummaryGapFilled{iModel} = titleStr;
            end

        end % iMethod
    end % iModel

    % Plot the annual mean data for all methods using the custom plotting function
    plotOceanVariableMaps(nppModelledAnnualMeanAll,qLon,qLat,myColourMap,cbString,...
        caxisMin,caxisMax,isCommonColourBar,globalNppStockSummary,'npp_annual_comp',[])

    % Plot the annual mean data for gap-filled data using the custom plotting function
    plotOceanVariableMaps(nppModelledAnnualMeanGapFilled,qLon,qLat,myColourMap,cbString,...
        caxisMin,caxisMax,isCommonColourBar,globalNppStockSummaryGapFilled,'npp_annual_gapfilled',[])
