
% ======================================================================= %
%                                                                         %
%            PAR0 climatology from NASA Aqua-MODIS and SeaWiFS            % 
%                                                                         %
% This script reads in photosynthetic available radiation at the surface  %
% ocean (PAR0) from NASA's Ocean Color website (Aqua-MODIS and SeaWiFS    %
% sensors).                                                               %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://oceancolor.gsfc.nasa.gov/l3/order/                          % 
%   - Time resolution: L3 monthly climatology                             %
%   - Space resolution: 4 km (8640 x 4320) (Aqua-MODIS),                  %                
%                       9 km (4320 x 2160) (SeaWiFS)                      %
%   - Version: 2002–2024 (Aqua-MODIS), 1997-2010 (SeaWiFS)                %
%   - Units: einstein m-2 d-1 = mol phot. m-2 d-1 (input), W m-2 (output) %
%                                                                         %
% Notice that if you click on the dataset link above, for the Aqua-MODIS  %
% sensor mapped monthly climatology, two months (July & August) are       % 
% missing from the collection of files. I accessed them following these   %
% instructions: https://forum.earthdata.nasa.gov/viewtopic.php?t=4784.    %
%                                                                         %
% This script uses this external function:                                %
%       processNASAsensorData.m - custom function                         %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 4 Nov 2024                                    %
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

fullpathInputPar0AquaModisDir = fullfile('data','raw','PAR0','AquaMODIS');
fullpathInputPar0SeawifsDir = fullfile('data','raw','PAR0','SeaWiFS');
fullpathOutputPar0AquaModisFile = fullfile('data','processed','par0_aquamodis.mat');
fullpathOutputPar0SeawifsFile = fullfile('data','processed','par0_seawifs.mat');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA FOR AQUA-MODIS
% -------------------------------------------------------------------------

S = ncinfo(fullfile(fullpathInputPar0AquaModisDir,... 
    ['AQUA_MODIS.' '20020901_20220930.L3m.MC.PAR.par.4km' '.nc']));

% Process data
[par0,lat,lon] = processSensorDataFromNASA(fullpathInputPar0AquaModisDir,...
    'AQUA_MODIS.','par','20020901_20220930.L3m.MC.PAR.par.4km');

% Units conversion: mol photons m-2 d-1 -> umol photons m-2 s-1 -> W m-2
par0_Wm2 = par0 .* (1e6 ./ (3600 .* 24)) .* (3.90e-19 .* 6.02e23 ./ 1e6); 

% Check for spurious data points
figure(); histogram(par0_Wm2(:), 100);

% Save the data
par0_lon = lon;
par0_lat = lat;
par0 = par0_Wm2;
save(fullpathOutputPar0AquaModisFile,'par0','par0_lat','par0_lon','-v7.3');

% Visual inspection
prepareDataForPlotting(fullpathOutputPar0AquaModisFile,[],'W m^{-2}',...
    0,200,true,'fig_monthly_par0_aquamodis','PAR0 Aqua-MODIS')
 
% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - EXTRACT THE DATA FOR SEAWIFS
% -------------------------------------------------------------------------

S = ncinfo(fullfile(fullpathInputPar0SeawifsDir,...
    ['SEASTAR_SEAWIFS_GAC.' '19980301_20100331.L3m.MC.PAR.par.9km' '.nc']));

% Process data
[par0,lat,lon] = processSensorDataFromNASA(fullpathInputPar0SeawifsDir,...
    'SEASTAR_SEAWIFS_GAC.','par','19980301_20100331.L3m.MC.PAR.par.9km');

% Units conversion: mol photons m-2 d-1 -> umol photons m-2 s-1 -> W m-2
par0_Wm2 = par0 .* (1e6 ./ (3600 .* 24)) .* (3.90e-19 .* 6.02e23 ./ 1e6); 

% Check for spurious data points
figure(); histogram(par0_Wm2(:), 100);

% Save the data
par0_lon = lon;
par0_lat = lat;
par0 = par0_Wm2;
save(fullpathOutputPar0SeawifsFile,'par0','par0_lat','par0_lon','-v7.3');

% Visual inspection
prepareDataForPlotting(fullpathOutputPar0SeawifsFile,[],'W m^{-2}',...
    0,200,true,'fig_monthly_par0_seawifs','PAR0 SeaWiFS')
 