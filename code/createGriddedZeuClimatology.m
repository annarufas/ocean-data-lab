
% ======================================================================= %
%                                                                         %
%                   Zeu climatology calculated from kd                    %
%                     from CMEMS and NASA Aqua-MODIS                      %
%                           and MLD from CMEMS                            % 
%                                                                         %
% This script creates a global gridded climatology of euphotic layer      %
% depth (zeu) using the diffuse attenuation coefficient at 490 nm (kd)    %
% product from Copernicus Marine Service (CMEMS) or NASA's Aqua-MODIS and %
% the mixed layer depth product from CMEMS. The output array is           %
% 1080 x 2160 x 12 and has units of m.                                    %                                        
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
%   Version 1.1 - 14 Nov 2024 (added NASA Aqua-MODIS)                     %
%                                                                         %
% ======================================================================= %

% Clear workspace, close figures, and add paths to plotting resources
close all; clear all; clc
addpath(genpath(fullfile('code')))
addpath(genpath(fullfile('resources','external'))); 
addpath(genpath(fullfile('resources','internal')));
addpath(genpath(fullfile('figures')))

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Choice of PAR0 level
isOnePercentPar0 = 0; % 1=zeu calculated as 1% PAR0, 0=zeu calculated as 0.1% PAR0

% Output filename
if isOnePercentPar0
    filenameTag = 'onepercentpar0';
else
    filenameTag = 'pointonepercentpar0';
end
fullpathOutputDir = fullfile('data','processed');

% Input datasets
fullpathInputKdCmems = fullfile(fullpathOutputDir,'kd_cmems_bgc.mat');
fullpathInputKdAquaModis = fullfile(fullpathOutputDir,'kd_modis.mat');
fullpathInputMldCmems = fullfile(fullpathOutputDir,'mld_cmems_phys.mat');

% Query points for interpolation (follow the BGC CMEMS grid, the one with the
% lowest resolution of the above)
qLats = linspace(-90, 90, 1080)';
qLons = linspace(-180, 180, 2160)';
   
% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD INPUT DATA AND REGRID TO COMMON GRID
% -------------------------------------------------------------------------

% Load input data
cmems = load(fullpathInputKdCmems,'kd','kd_lat','kd_lon');
aquamodis = load(fullpathInputKdAquaModis,'kd','kd_lat','kd_lon') ;
load(fullpathInputMldCmems,'mld','mld_lat','mld_lon');

% Data grid
[Xc, Yc, Tc] = ndgrid(cmems.kd_lat, cmems.kd_lon, (1:12)');
[Xa, Ya, Ta] = ndgrid(aquamodis.kd_lat, aquamodis.kd_lon, (1:12)');
[Xm, Ym, Tm] = ndgrid(mld_lat, mld_lon, (1:12)');

% Query grid
[qX, qY, qT] = ndgrid(qLats, qLons, (1:12)');

% Interpolant -use first-order (linear) interpolation and extrapolation
Fkdc = griddedInterpolant(Xc, Yc, Tc, cmems.kd, 'linear'); 
Fkda = griddedInterpolant(Xa, Ya, Ta, aquamodis.kd, 'linear');
Fmld = griddedInterpolant(Xm, Ym, Tm, mld, 'linear'); 

% Regrid 
qKdCmems = Fkdc(qX, qY, qT);
qKdAquaModis = Fkda(qX, qY, qT);
qMld = Fmld(qX, qY, qT);
% figure(); pcolor(qKdAquaModis(:,:,4)); caxis([0 0.2]); shading interp;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE ZEU USING CMEMS KD AND MLD FROM CMEMS
% -------------------------------------------------------------------------

[zeu,~] = calculateZeuFromKdAndMLD(qKdCmems,qMld,isOnePercentPar0);
                 
% Check for spurious data points
%figure(); histogram(zeu(:), 100);

% Save output
zeu_lat = qLats;
zeu_lon = qLons;
fullpathOutputZeuFile = fullfile(fullpathOutputDir,strcat('zeu_calculated_kdcmems_mldcmems_',filenameTag,'.mat'));
save(fullfile(fullpathOutputZeuFile),'zeu','zeu','zeu_lat','zeu_lon','-v7.3')   

% Visual inspection
prepareDataForPlotting(fullpathOutputZeuFile,[],'m',...
    0,200,true,strcat('fig_monthly_zeu_calculated_kdcmems_mldcmems_',filenameTag),'Zeu calculated from CMEMS kd')

clear zeu

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CALCULATE ZEU USING AQUA-MODIS KD AND MLD FROM CMEMS
% -------------------------------------------------------------------------

[zeu,~] = calculateZeuFromKdAndMLD(qKdAquaModis,qMld,isOnePercentPar0);        
                  
% Check for spurious data points
%figure(); histogram(zeu(:), 100);

% Save output
zeu_lat = qLats;
zeu_lon = qLons;
fullpathOutputZeuFile = fullfile(fullpathOutputDir,strcat('zeu_calculated_kdmodis_mldcmems_',filenameTag,'.mat'));
save(fullfile(fullpathOutputZeuFile),'zeu','zeu','zeu_lat','zeu_lon','-v7.3')   

% Visual inspection
prepareDataForPlotting(fullpathOutputZeuFile,[],'m',...
    0,200,true,strcat('fig_monthly_zeu_calculated_kdmodis_mldcmems_',filenameTag),'Zeu calculated from Aqua-MODIS kd')
