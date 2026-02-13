
% ======================================================================= %
%                                                                         %
%                   Zeu climatology calculated from kd                    %
%                     from CMEMS or NASA Aqua-MODIS /                     %  
%                           chla from OC-CCI                              %
%                    and MLD from CMEMS or IFREMER                        % 
%                                                                         %
% This script creates a global gridded climatology of euphotic layer      %
% depth (zeu) using the mixed layer depth product from CMEMS or IFREMER,  %
% and either (i) diffuse attenuation coefficient at 490 nm (kd) product   %
% from Copernicus Marine Service (CMEMS) or NASA's Aqua-MODIS, or (ii)    %
% chlorphyll concentration product from OC-CCI. The output array is       %
% 1080 x 2160 x 12 and has units of m.                                    %                                        
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
%   Version 1.1 - 14 Nov 2024 (added NASA Aqua-MODIS)                     %
%   Version 1.2 - 23 Apr 2024 (added Chla OC-CCI and MLD IFREMER)         %
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
fullpathInputChlaOccci = fullfile(fullpathOutputDir,'chla_occci.mat');
fullpathInputMldCmems = fullfile(fullpathOutputDir,'mld_cmems_phys.mat');
fullpathInputMldIfremer = fullfile(fullpathOutputDir,'mld_ifremer.mat');

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
kdcmems = load(fullpathInputKdCmems,'kd','kd_lat','kd_lon');
kdaquamodis = load(fullpathInputKdAquaModis,'kd','kd_lat','kd_lon') ;
load(fullpathInputChlaOccci,'chla','chla_lat','chla_lon');
mldcmems = load(fullpathInputMldCmems,'mld','mld_lat','mld_lon');
mldifremer = load(fullpathInputMldIfremer,'mld','mld_lat','mld_lon');

% Data grid
[Xc, Yc, Tc] = ndgrid(kdcmems.kd_lat, kdcmems.kd_lon, (1:12)');
[Xa, Ya, Ta] = ndgrid(kdaquamodis.kd_lat, kdaquamodis.kd_lon, (1:12)');
[Xcl, Ycl, Tcl] = ndgrid(chla_lat, chla_lon, (1:12)');
[Xmc, Ymc, Tmc] = ndgrid(mldcmems.mld_lat, mldcmems.mld_lon, (1:12)');
[Xmi, Ymi, Tmi] = ndgrid(mldifremer.mld_lat, mldifremer.mld_lon, (1:12)');

% Query grid
[qX, qY, qT] = ndgrid(qLats, qLons, (1:12)');

% Interpolant -use first-order (linear) interpolation and extrapolation
Fkdc = griddedInterpolant(Xc, Yc, Tc, kdcmems.kd, 'linear'); 
Fkda = griddedInterpolant(Xa, Ya, Ta, kdaquamodis.kd, 'linear');
Fchla = griddedInterpolant(Xcl, Ycl, Tcl, chla, 'linear');
Fmldc = griddedInterpolant(Xmc, Ymc, Tmc, mldcmems.mld, 'linear'); 
Fmldi = griddedInterpolant(Xmi, Ymi, Tmi, mldifremer.mld, 'linear'); 

% Regrid 
qKdCmems = Fkdc(qX, qY, qT);
qKdAquaModis = Fkda(qX, qY, qT);
qChla = Fchla(qX, qY, qT);
qMldCmems = Fmldc(qX, qY, qT);
qMldIfremer = Fmldi(qX, qY, qT);
% figure(); pcolor(qKdAquaModis(:,:,4)); caxis([0 0.2]); shading interp;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE ZEU USING CMEMS KD AND MLD FROM CMEMS
% -------------------------------------------------------------------------

[zeu,~] = calculateZeuFromKdAndMLD(qKdCmems,qMldCmems,isOnePercentPar0);
             
% Check for spurious data points
%figure(); histogram(zeu(:), 100);

% Save output
zeu_lat = qLats;
zeu_lon = qLons;
fullpathOutputZeuFile = fullfile(fullpathOutputDir,strcat('zeu_calculated_kdcmems_mldcmems_',filenameTag,'.mat'));
save(fullfile(fullpathOutputZeuFile),'zeu','zeu_lat','zeu_lon','-v7.3')   

% Visual inspection
prepareDataForPlotting(fullpathOutputZeuFile,[],'m',...
    0,200,true,strcat('fig_monthly_zeu_calculated_kdcmems_mldcmems_',filenameTag),'Zeu calculated from CMEMS kd')

clear zeu

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CALCULATE ZEU USING AQUA-MODIS KD AND MLD FROM CMEMS
% -------------------------------------------------------------------------

[zeu,~] = calculateZeuFromKdAndMLD(qKdAquaModis,qMldCmems,isOnePercentPar0);        
                  
% Check for spurious data points
%figure(); histogram(zeu(:), 100);

% Save output
zeu_lat = qLats;
zeu_lon = qLons;
fullpathOutputZeuFile = fullfile(fullpathOutputDir,strcat('zeu_calculated_kdmodis_mldcmems_',filenameTag,'.mat'));
save(fullfile(fullpathOutputZeuFile),'zeu','zeu_lat','zeu_lon','-v7.3')   

% Visual inspection
prepareDataForPlotting(fullpathOutputZeuFile,[],'m',...
    0,200,true,strcat('fig_monthly_zeu_calculated_kdmodis_mldcmems_',filenameTag),'Zeu calculated from Aqua-MODIS kd')

clear zeu

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - CALCULATE ZEU USING OC-CCI CHLA AND MLD FROM IFREMER
% -------------------------------------------------------------------------

[zeu,~] = calculateZeuFromChlaAndMLD(qChla,qMldIfremer,isOnePercentPar0);        
                  
% Check for spurious data points
%figure(); histogram(zeu(:), 100);

% Save output
zeu_lat = qLats;
zeu_lon = qLons;
fullpathOutputZeuFile = fullfile(fullpathOutputDir,strcat('zeu_calculated_chlaoccci_mldifremer_',filenameTag,'.mat'));
save(fullfile(fullpathOutputZeuFile),'zeu','zeu_lat','zeu_lon','-v7.3')   

% Visual inspection
prepareDataForPlotting(fullpathOutputZeuFile,[],'m',...
    0,200,true,strcat('fig_monthly_zeu_calculated_chlaoccci_mldifremer_',filenameTag),'Zeu calculated from OC-CCI Chla')

clear zeu