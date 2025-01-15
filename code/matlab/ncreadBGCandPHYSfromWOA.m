
% ======================================================================= %
%                                                                         %
%     Nitrate, phosphate, silicate, oxygen, temperature and salinity      %
%                       climatologies from WOA                            %
%                                                                         %
% This script reads in World Ocean Atlas 2023 (WOA23) monthly and annual  %
% climatologies of BGC (nitrate, phosphate, silicate and oxygen) and PHYS % 
% (temperature, salinity) variables and merges the two at depth to create %
% a comprehensive monthly climatology.                                    %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/            %
%   - Time resolution: monthly and annual climatologies                   %
%   - Space resolution: 1° x 1° (180 x 360)                               %
%   - Version: 2023                                                       %
%   - Units: Nitrate: mmol m-3 (=umol kg-1)                               % 
%            Phosphate: mmol m-3                                          %
%            Silicate: mmol m-3                                           % 
%            Oxygen: mL L-1 (=umol kg-1)                                  %
%            Temperature: °C                                              % 
%            Salinity: PSU                                                %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 15 Nov 2024                                   %
%   Version 1.1 - 3 Jan 2025: added salinity and density calculation      %
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

fullpathInputWoaDir = fullfile('data','raw','WOA','WOA23');
fullpathOutputWoaDir = fullfile('data','processed');

% WOA version
commonTag = 'woa23'; % for World Ocean Atlas 2023 version

% Inspection for variable names
S1 = ncinfo(fullfile(fullpathInputWoaDir, [commonTag '_decav_t01_01.nc'])); 
S2 = ncinfo(fullfile(fullpathInputWoaDir, [commonTag '_all_n01_01.nc']));

% Datasets metadata {variable name, file infix, NetCDF variable name}
datasetsMetadata = {...
    'nit',  'all_n',   'n_an';...
    'phos', 'all_p',   'p_an';...
    'sil',  'all_i',   'i_an';...
    'oxy',  'all_o',   'o_an';...
    'temp', 'decav_t', 't_an'; ...
    'sal',  'decav_s', 's_an'};

% Ouput files that will be created by this script:

% Monthly climatologies
filenameOutputWoaMonthlyNit = ['nit_monthly_' commonTag '.mat'];
filenameOutputWoaMonthlyPhos = ['phos_monthly_' commonTag '.mat'];
filenameOutputWoaMonthlySil = ['sil_monthly_' commonTag '.mat'];
filenameOutputWoaMonthlyOxy = ['oxy_monthly_' commonTag '.mat'];
filenameOutputWoaMonthlyTemp = ['temp_monthly_' commonTag '.mat'];
filenameOutputWoaMonthlySal = ['sal_monthly_' commonTag '.mat'];

% Annual climatology
filenameOutputWoaAnnualTemp = ['temp_annual_' commonTag '.mat'];

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - GET GRID DIMENSIONS NAD VARIABLE NAMING
% -------------------------------------------------------------------------

% Read longitude and latitude (common to all variables and monthly and annual data)
ncFile = fullfile(fullpathInputWoaDir, [commonTag '_all_o01_01.nc']);
lon = double(ncread(ncFile, 'lon'));
lat = double(ncread(ncFile, 'lat'));

% Read metadata
varNames = datasetsMetadata(:,1);
ncFileInfix = datasetsMetadata(:,2);
ncFileVarName = datasetsMetadata(:,3);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - EXTRACT MONTHLY DATA
% -------------------------------------------------------------------------

% Read depth levels for each variable (first month, 01, same for all)
suffix = '01_01';
[~,monthlyDataDepths] = readWoaDepthLevels(varNames,ncFileInfix,fullpathInputWoaDir,commonTag,suffix);

% Read monthly data
monthlyData = struct();
for iMonth = 1:12
    for iVar = 1:numel(varNames)
        ncFile = fullfile(fullpathInputWoaDir,[commonTag '_' ncFileInfix{iVar} sprintf('%02d',iMonth) '_01.nc']);
        monthlyData.(varNames{iVar})(:,:,:,iMonth) = double(ncread(ncFile, ncFileVarName{iVar}));
    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - EXTRACT ANNUAL DATA
% -------------------------------------------------------------------------

% Read depth levels for each variable (first month, 01, same for all)
suffix = '00_01';
[~,annualDataDepths] = readWoaDepthLevels(varNames,ncFileInfix,fullpathInputWoaDir,commonTag,suffix);

% Read annual data
annualData = struct();
for iVar = 1:numel(varNames)
    ncFile = fullfile(fullpathInputWoaDir,[commonTag '_' ncFileInfix{iVar} '00_01.nc']);
    annualData.(varNames{iVar}) = double(ncread(ncFile, ncFileVarName{iVar}));
end

% Read standard deviation for annual temperature
annualTempStd = double(ncread(fullfile(fullpathInputWoaDir, [commonTag '_decav_t00_01.nc']),'t_sd'));

% % Visual inspection
% iNitAnnual100m = find(annualDataDepths.nit == 100);
% iPhosAnnual100m = find(annualDataDepths.phos == 100);
% iOxyAnnual100m = find(annualDataDepths.oxy == 100);
% iTempAnnual100m = find(annualDataDepths.temp == 100);
% 
% figure(); pcolor(flipud(rot90(annualData.temp(:,:,iTempAnnual100m)))); 
% caxis([-2 25]); shading interp; colormap(jet);
% 
% figure(); pcolor(flipud(rot90(annualData.nit(:,:,iNitAnnual100m)))); 
% caxis([0 27]); shading interp; colormap(jet)
% 
% figure(); pcolor(flipud(rot90(annualData.oxy(:,:,iOxyAnnual100m)))); 
% caxis([0 7]); shading interp; colormap(jet)
% 
% figure(); pcolor(flipud(rot90(annualData.phos(:,:,iPhosAnnual100m)))); 
% caxis([0 7]); shading interp; colormap(jet)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - ADD DEEP ANNUAL DATA TO MONTHLY ARRAYS
% -------------------------------------------------------------------------

dataCombined = struct();
dataCombinedDepths = struct();

for iVar = 1:numel(varNames)
    varName = varNames{iVar};
    
    % Find last depth populated in monthly array
    iLastDepthMonthly = find(annualDataDepths.(varName) == monthlyDataDepths.(varName)(end));

    % Extract the deep annual data for this variable beyond the last monthly depth
    deepAnnualData = annualData.(varName)(:,:,(iLastDepthMonthly+1):end);
    
    % Concatenate the annual and monthly depth profiles
    for iMonth = 1:12
        dataCombined.(varName)(:,:,:,iMonth) = ...
            cat(3, monthlyData.(varName)(:,:,:,iMonth), deepAnnualData);
    end
    
    % Get the depths
     dataCombinedDepths.(varName) = ...
        cat(1, monthlyDataDepths.(varName), annualDataDepths.(varName)((iLastDepthMonthly+1):end));
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - CHECK AND SAVE THE DATA
% -------------------------------------------------------------------------

% Check for spurious data points
% figure(); histogram(dataCombined.nit(:),100);
% figure(); histogram(dataCombined.phos(:),100);
% figure(); histogram(dataCombined.sil(:),100);
% figure(); histogram(dataCombined.oxy(:),100);
% figure(); histogram(dataCombined.temp(:),100);
% figure(); histogram(dataCombined.sal(:),100);

% Permute data (from lon x lat to lat x lon)
permutedData = struct();
for iVar = 1:numel(varNames)
    permutedData.(varNames{iVar}) = permute(dataCombined.(varNames{iVar}), [2, 1, 3, 4]);
end
permutedData.temp_annual = permute(annualData.temp, [2, 1, 3, 4]);
permutedData.temp_annual_std = permute(annualTempStd, [2, 1, 3, 4]);

% Output arrays naming
woa_lat = lat;
woa_lon = lon;
woa_annual_lat = lat;
woa_annual_lon = lon;

nit = permutedData.nit;
phos = permutedData.phos;
sil = permutedData.sil;
oxy = permutedData.oxy;
temp = permutedData.temp;
sal = permutedData.sal;
temp_annual = permutedData.temp_annual;
temp_annual_std = permutedData.temp_annual_std;

woa_depth_nit  = dataCombinedDepths.nit;
woa_depth_phos = dataCombinedDepths.phos;
woa_depth_sil  = dataCombinedDepths.sil;
woa_depth_oxy  = dataCombinedDepths.oxy;
woa_depth_temp = dataCombinedDepths.temp;
woa_depth_sal  = dataCombinedDepths.sal;
woa_annual_depth_temp = annualDataDepths.temp;

% % Visual evolution over depth
% figure()
% myVideo = VideoWriter(fullfile(fullpathOutputWoaDir,'NitEvolThroughDepth.avi'));
% myVideo.FrameRate = 5;  
% myVideo.Quality = 100;    
% open(myVideo)
% for iDepth = 1:numel(woa_depth_nit)
%     pcolor(flipud(rot90(nit(:,:,iDepth,5)))); 
%     caxis([0 40]);
%     cb = colorbar('FontSize', 12); 
%     cb.Label.String = 'NO3, mmol m-3';
%     shading interp
%     colormap(jet)
%     box on
%     frame = getframe(gcf);
%     writeVideo(myVideo,frame); 
% end
% close(myVideo)

% Save combined monthly + annual arrays
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyNit),...
    'nit','woa_lat','woa_lon','woa_depth_nit')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyPhos),...
    'phos','woa_lat','woa_lon','woa_depth_phos')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlySil),...
    'sil','woa_lat','woa_lon','woa_depth_sil')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyOxy),...
    'oxy','woa_lat','woa_lon','woa_depth_oxy')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyTemp),...
    'temp','woa_lat','woa_lon','woa_depth_temp')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlySal),...
    'sal','woa_lat','woa_lon','woa_depth_sal')

% Save annual array
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaAnnualTemp),...
    'temp_annual','temp_annual_std','woa_annual_depth_temp','woa_annual_lat','woa_annual_lon')

% Visual inspection
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyNit),11,'mmol m^{-3}',...
    0,25,true,'fig_monthly_nit_woa','Nitrate at 50 m, WOA23')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyPhos),11,'mmol m^{-3}',...
    0,2.5,true,'fig_monthly_phos_woa','Phosphate at 50 m, WOA23')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlySil),11,'mmol m^{-3}',...
    0,50,true,'fig_monthly_sil_woa','Silicate at 50 m, WOA23')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyOxy),11,'mL L^{-1}',...
    160,360,true,'fig_monthly_oxy_woa','Oxygen at 50 m, WOA23')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyTemp),11,'°C',...
    0,25,true,'fig_monthly_temp_woa','Temperature at 50 m, WOA23')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlySal),11,'PSU',...
    32,36,true,'fig_monthly_sal_woa','Salinity at 50 m, WOA23')

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function [nz,depths] = readWoaDepthLevels(varNames,ncFileInfix,fullpathInputWoaDir,commonTag,suffix)

    depths = cell2struct(cell(1, numel(varNames)), varNames, 2);

    % Populate the depths structure
    for iVar = 1:numel(varNames)
        varName = varNames{iVar};
        ncFile = fullfile(fullpathInputWoaDir, [commonTag '_' ncFileInfix{iVar} suffix '.nc']);
        depths.(varName) = double(ncread(ncFile, 'depth'));
    end

    % Get the number of depth levels for each variable
    nz = structfun(@length, depths);

end % readWoaDepthLevels
