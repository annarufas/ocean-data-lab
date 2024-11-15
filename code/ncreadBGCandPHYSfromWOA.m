
% ======================================================================= %
%                                                                         %
%         Nitrate, phosphate, silicate, oxygen and temperature            %
%                       climatologies from WOA                            %
%                                                                         %
% This script reads in World Ocean Atlas 2023 (WOA23) monthly and annual  %
% climatologies of BGC (nitrate, phosphate, silicate and oxygen) and PHYS % 
% (temperature) variables and merges the two at depth to create a         %
% comprehensive monthly climatology.                                      %
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
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 15 Nov 2024                                   %
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

filenameInputCommonTag = 'woa23_';
fullpathInputWoaDir = fullfile('data','raw','WOA','WOA23');
fullpathOutputWoaDir = fullfile('data','processed');
filenameOutputWoaMonthlyNit = 'nit_monthly_woa23.mat';
filenameOutputWoaMonthlyPhos = 'phos_monthly_woa23.mat';
filenameOutputWoaMonthlySil = 'sil_monthly_woa23.mat';
filenameOutputWoaMonthlyTemp = 'temp_monthly_woa23.mat';
filenameOutputWoaMonthlyOxy = 'oxy_monthly_woa23.mat';
filenameOutputWoaAnnualTemp = 'temp_annual_woa23.mat';

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

S1 = ncinfo(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'decav_t01_01.nc'])); % monthly
S2 = ncinfo(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_n01_01.nc']));
S3 = ncinfo(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'decav_t00_01.nc'])); % annual
S4 = ncinfo(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_n00_01.nc'])); 

% Monthly data
lon = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_o01_01.nc']),'lon')); 
lat = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_o01_01.nc']),'lat'));
zTemp = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'decav_t01_01.nc']),'depth'));
zNit = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_n01_01.nc']),'depth'));
zOxy = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_o01_01.nc']),'depth'));
zSil = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_i01_01.nc']),'depth'));
zPhos = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_p01_01.nc']),'depth'));

% Get sizes of data arrays
nx = length(lon);
ny = length(lat);
nzTemp = length(zTemp);
nzNit = length(zNit);
nzOxy = length(zOxy);
nzSil = length(zSil);
nzPhos = length(zPhos);

% Initialise monthly arrays
monthlyTemp = zeros([nx ny nzTemp 12]);
monthlyNit = zeros([nx ny nzNit 12]);
monthlyOxy = zeros([nx ny nzOxy 12]);
monthlySil = zeros([nx ny nzSil 12]);
monthlyPhos = zeros([nx ny nzPhos 12]);

for iMonth = 1:12
   monthlyTemp(:,:,:,iMonth) = ncread(fullfile(fullpathInputWoaDir,...
       [filenameInputCommonTag 'decav_t' sprintf('%02d',iMonth) '_01.nc']), 't_an');
   monthlyNit(:,:,:,iMonth) = ncread(fullfile(fullpathInputWoaDir,...
       [filenameInputCommonTag 'all_n' sprintf('%02d',iMonth) '_01.nc']), 'n_an');
   monthlyOxy(:,:,:,iMonth) = ncread(fullfile(fullpathInputWoaDir,...
       [filenameInputCommonTag 'all_o' sprintf('%02d',iMonth) '_01.nc']), 'o_an');   
   monthlySil(:,:,:,iMonth) = ncread(fullfile(fullpathInputWoaDir,...
       [filenameInputCommonTag 'all_i' sprintf('%02d',iMonth) '_01.nc']), 'i_an');
   monthlyPhos(:,:,:,iMonth) = ncread(fullfile(fullpathInputWoaDir,...
       [filenameInputCommonTag 'all_p' sprintf('%02d',iMonth) '_01.nc']), 'p_an');
end

% .........................................................................

% Annual data

lonA = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_o00_01.nc']),'lon'));
latA = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_o00_01.nc']),'lat'));
zTempA = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'decav_t00_01.nc']),'depth'));
zNitA = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_n00_01.nc']),'depth'));
zOxyA = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_o00_01.nc']),'depth'));
zSilA = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_i00_01.nc']),'depth'));
zPhosA = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_p00_01.nc']),'depth'));

annualTemp = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'decav_t00_01.nc']),'t_an'));
annualNit = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_n00_01.nc']),'n_an'));
annualOxy = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_o00_01.nc']),'o_an'));
annualSil = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_i00_01.nc']),'i_an'));
annualPhos = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'all_p00_01.nc']),'p_an'));

% Read standard deviation for annual temperature
annualTempStd = double(ncread(fullfile(fullpathInputWoaDir, [filenameInputCommonTag 'decav_t00_01.nc']),'t_sd'));

% % Visual inspection
% iTempAnnual100m = find(zTempA == 100);
% iNitAnnual100m = find(zNitA == 100);
% iOxyAnnual100m = find(zOxyA == 100);
% iPhosAnnual100m = find(zPhosA == 100);
% 
% figure()
% pcolor(flipud(rot90(annualTemp(:,:,iTempAnnual100m)))); 
% caxis([-2 25]); 
% cb = colorbar('FontSize', 15, 'FontWeight', 'bold'); 
% cb.Label.String = 'Temperature (°C)';
% shading interp
% colormap(jet)
% box on
% 
% figure()
% pcolor(flipud(rot90(annualNit(:,:,iNitAnnual100m)))); 
% caxis([0 27]); 
% cb = colorbar('FontSize', 15, 'FontWeight', 'bold'); 
% cb.Label.String = 'Nitrate (mmol m^{-3})';
% shading interp
% colormap(jet)
% box on
% 
% figure()
% pcolor(flipud(rot90(annualOxy(:,:,iOxyAnnual100m)))); 
% % caxis([0 7]); 
% cb = colorbar('FontSize', 15, 'FontWeight', 'bold'); 
% cb.Label.String = 'Oxygen (ml L^{-1})';
% shading interp
% colormap(jet)
% box on
% 
% figure()
% pcolor(flipud(rot90(annualPhos(:,:,iPhosAnnual100m)))); 
% % caxis([0 7]); 
% cb = colorbar('FontSize', 15, 'FontWeight', 'bold'); 
% cb.Label.String = 'Phosphate (mmol m^{-3})';
% shading interp
% colormap(jet)
% box on

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - ADD DEEP ANNUAL CLIMATOLOGICAL DATA TO MONTHLY CLIMATOLOGICAL
% ARRAYS
% -------------------------------------------------------------------------

% Crop data from the annual mean
% Monthly oxygen and temperature data go down to depth of 1500 m
% Monthly nitrate, phosphate and silicate data go down to depth of 800 m
iLastDepthMonthlyTemp = find(zTempA == zTemp(end));
iLastDepthMonthlyNit  = find(zNitA == zNit(end));
iLastDepthMonthlyOxy  = find(zOxyA == zOxy(end));
iLastDepthMonthlySil  = find(zSilA == zSil(end));
iLastDepthMonthlyPhos = find(zPhosA == zPhos(end));

deepAnnualTemp = annualTemp(:,:,(iLastDepthMonthlyTemp+1:end));
deepAnnualNit  = annualNit(:,:,(iLastDepthMonthlyNit+1:end));
deepAnnualOxy  = annualOxy(:,:,(iLastDepthMonthlyOxy+1:end)); 
deepAnnualSil  = annualSil(:,:,(iLastDepthMonthlySil+1:end)); 
deepAnnualPhos = annualPhos(:,:,(iLastDepthMonthlyPhos+1:end)); 

% Concatenate the matrices
nzTempWaterCol = nzTemp + length(deepAnnualTemp(1,1,:));
nzNitWaterCol  = nzNit + length(deepAnnualNit(1,1,:));
nzOxyWaterCol  = nzOxy + length(deepAnnualOxy(1,1,:));
nzSilWaterCol  = nzSil + length(deepAnnualSil(1,1,:));
nzPhosWaterCol = nzPhos + length(deepAnnualPhos(1,1,:));

allTemp = zeros([nx ny nzTempWaterCol 12]);
allNit = zeros([nx ny nzNitWaterCol 12]);
allOxy = zeros([nx ny nzOxyWaterCol 12]);
allSil = zeros([nx ny nzSilWaterCol 12]);
allPhos = zeros([nx ny nzPhosWaterCol 12]);

for iMonth = 1:12
    allTemp(:,:,:,iMonth) = cat(3,monthlyTemp(:,:,:,iMonth),deepAnnualTemp(:,:,:));
    allNit(:,:,:,iMonth) = cat(3,monthlyNit(:,:,:,iMonth),deepAnnualNit(:,:,:));
    allOxy(:,:,:,iMonth) = cat(3,monthlyOxy(:,:,:,iMonth),deepAnnualOxy(:,:,:));
    allSil(:,:,:,iMonth) = cat(3,monthlySil(:,:,:,iMonth),deepAnnualSil(:,:,:));
    allPhos(:,:,:,iMonth) = cat(3,monthlyPhos(:,:,:,iMonth),deepAnnualPhos(:,:,:));    
end

woa_depth_temp = cat(1,zTemp,zTempA(iLastDepthMonthlyTemp+1:end));
woa_depth_nit  = cat(1,zNit,zNitA(iLastDepthMonthlyNit+1:end));
woa_depth_oxy  = cat(1,zOxy,zOxyA(iLastDepthMonthlyOxy+1:end));
woa_depth_sil  = cat(1,zSil,zSilA(iLastDepthMonthlySil+1:end));
woa_depth_phos = cat(1,zPhos,zPhosA(iLastDepthMonthlyPhos+1:end));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CHECK AND SAVE THE DATA
% -------------------------------------------------------------------------

% Check for spurious data points
figure(); histogram(allTemp(:),100);
figure(); histogram(allNit(:),100);
figure(); histogram(allOxy(:),100);
figure(); histogram(allSil(:),100);
figure(); histogram(allPhos(:),100);

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

% Permute dimensions arrangement before saving
nit = permute(allNit, [2, 1, 3, 4]); 
phos = permute(allPhos, [2, 1, 3, 4]); 
sil = permute(allSil, [2, 1, 3, 4]); 
temp = permute(allTemp, [2, 1, 3, 4]); 
oxy = permute(allOxy, [2, 1, 3, 4]); 
woa_lat = lat;
woa_lon = lon;
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyNit),...
    'nit','woa_lat','woa_lon','woa_depth_nit')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyPhos),...
    'phos','woa_lat','woa_lon','woa_depth_phos')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlySil),...
    'sil','woa_lat','woa_lon','woa_depth_sil')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyTemp),...
    'temp','woa_lat','woa_lon','woa_depth_temp')
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyOxy),...
    'oxy','woa_lat','woa_lon','woa_depth_oxy')

% Permute dimensions arrangement before saving
woa_annual_lat = latA;
woa_annual_lon = lonA;
woa_annual_depth_temp = zTempA;
temp_annual = permute(annualTemp, [2, 1, 3]);
temp_annual_std = permute(annualTempStd, [2, 1, 3]);
save(fullfile(fullpathOutputWoaDir,filenameOutputWoaAnnualTemp),...
    'temp_annual','temp_annual_std','woa_annual_depth_temp','woa_annual_lat','woa_annual_lon')

% Visual inspection
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyNit),11,'mmol m^{-3}',...
    0,25,true,'fig_monthly_nit_woa','Nitrate at 50 m, WOA')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyPhos),11,'mmol m^{-3}',...
    0,2.5,true,'fig_monthly_phos_woa','Phosphate at 50 m, WOA')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlySil),11,'mmol m^{-3}',...
    0,50,true,'fig_monthly_sil_woa','Silicate at 50 m, WOA')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyTemp),11,'°C',...
    0,25,true,'fig_monthly_temp_woa','Temperature at 50 m, WOA')
prepareDataForPlotting(fullfile(fullpathOutputWoaDir,filenameOutputWoaMonthlyOxy),11,'mL L^{-1}',...
    160,360,true,'fig_monthly_oxy_woa','Oxygen at 50 m, WOA')
