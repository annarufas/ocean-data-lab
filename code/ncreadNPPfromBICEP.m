
% ======================================================================= %
%                                                                         %
%                 NPP climatology from BICEP project                      %
%                                                                         %
% This script reads in net primary production (NPP) data from the         %
% Biological Pump and Carbon Exchange Processes (BICEP) project.          %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://catalogue.ceda.ac.uk/uuid/69b2c9c6c4714517ba10dab3515e4ee6/ % 
%   - Time resolution: monthly climatology                                %
%   - Space resolution: 9 km resolution (2160 x 4320)                     %
%   - Version: v4.2, 1998–2020                                            %
%   - Units: mg C m-2 d-1                                                 %
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

filenameOutputNpp = 'npp_bicep.mat';
fullpathInputNppDir = fullfile('data','raw','BICEP','BICEP_NPP');
fullpathOutputNppFile = fullfile('data','processed',filenameOutputNpp);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

% Time vector
yearsVector = 1998:2020;

for iMonth = 1:12

    for iYear = 1:numel(yearsVector)
        
        thisYear = yearsVector(iYear);
        yearFolderPath = fullfile(fullpathInputNppDir, num2str(thisYear));
        fileNames = dir(fullfile(yearFolderPath,'*.nc'));

        filePath = fullfile(yearFolderPath, fileNames(iMonth).name);
        %S = ncinfo(filePath); % short summary

        % Read longitude and latitude
        lat = ncread(filePath,'latitude');
        lon = ncread(filePath,'longitude');
        
        if (iYear == 1 && iMonth == 1)
            netprimaryproduction = NaN(numel(lat),numel(lon),12); 
            netprimaryproduction_err = NaN(numel(lat),numel(lon),12); 
        end

        % Read data and permute lat and lon
        D = zeros(numel(lat),numel(lon));
        Dtmp = ncread(filePath,'pp');
        Dperm = permute(Dtmp,[2 1]); % swap lat and lon
        D(:,:) = Dperm;

        % Sort lat and lon to have monotonically increasing values and 
        % apply the sorting to the data
        [lat_sort, sortLatIdx] = sort(lat);
        [lon_sort, sortLonIdx] = sort(lon);
        D_sort = D(sortLatIdx,sortLonIdx); % lat x lon x vars

        % Initialise the output data array on the first iteration
        if (iYear == 1)    
            D_out = NaN(numel(lat_sort),numel(lon_sort),numel(yearsVector),'single');
        end

        % Store the sorted data in the output array at the calculated time index
        D_out(:,:,iYear) = D_sort(:,:);

    end % iYear
    
    netprimaryproduction(:,:,iMonth) = mean(D_out(:,:,:),3,'omitnan'); % mg C m-2 d-1
    netprimaryproduction_err(:,:,iMonth) = std(D_out(:,:,:),0,3,'omitnan');
   
end % iMonth

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CHECK VALUES AND SAVE
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(netprimaryproduction(:),100);

% Save the data
npp_lon = lon_sort;
npp_lat = lat_sort;
npp_avg = netprimaryproduction;
npp_err = netprimaryproduction_err;
save(fullpathOutputNppFile,'npp_avg','npp_err','npp_lat','npp_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputNppFile,[],'mg C m-2 d-1',...
    0,1000,true,'fig_monthly_npp_bicep','NPP BICEP')
