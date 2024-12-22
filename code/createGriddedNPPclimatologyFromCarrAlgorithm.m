
% ======================================================================= %
%                                                                         %
%               NPP climatology from Carr 2002 algorithm                  % 
%                                                                         %
% This script creates a global gridded climatology of net primary         %
% production (NPP) using the algorithm of Carr 2002                       %
% (https://doi.org/10.1016/S0967-0645(01)00094-7) and input data of       %
% chlorophyll a (chla), sea surface temperature (SST) and photosynthetic  %
% active radiation in the surface ocean (PAR0). The output array is a     %
% monthly climatology of 180 x 360 x 12 pixels and has units of           %
% mg C m-2 d-1.                                                           %                                        
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 4 Nov 2024                                    %
%   Version 1.1 - 21 Dec 2024: added Aqua-MODIS products                  %
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

% Choice of euphotic layer depth for the implementation of Carr 2002
% algorithm, where 1=Carr 2002, 0=Behrenfeld & Falkowski 1997 (as in Henson et al. 2012)
isZeuCarr = 1; 

% Choice of PAR0 and chla from Aqua-MODIS or Seawifs sensor, 
% where 1=Aqua-MODIS, 0=Seawifs (as in Henson et al. 2012)
isSensorAquaModis = 1; 

% Input files
if isSensorAquaModis
    fullpathInputChla = fullfile('data','processed','chla_aquamodis.mat');  % mg m-3
    fullpathInputPar0 = fullfile('data','processed','par0_aquamodis.mat');  % W m-2
else ~isSensorAquaModis
    fullpathInputChla = fullfile('data','processed','chla_seawifs.mat');  % mg m-3
    fullpathInputPar0 = fullfile('data','processed','par0_seawifs.mat');  % W m-2
end
fullpathInputSst  = fullfile('data','processed','sst_pathfinder_v5.mat'); % ºC

% Output files
if isSensorAquaModis && isZeuCarr
    fullpathOutputNppFile = fullfile('data','processed','npp_carr2002_aquamodis_pathfinder_zeuc02.mat');
    fullpathFigureName = 'fig_monthly_npp_carr2002_aquamodis_zeuc02';
elseif isSensorAquaModis && ~isZeuCarr
    fullpathOutputNppFile = fullfile('data','processed','npp_carr2002_aquamodis_pathfinder_zeub97.mat');
    fullpathFigureName = 'fig_monthly_npp_carr2002_aquamodis_zeub97';
elseif ~isSensorAquaModis && isZeuCarr
    fullpathOutputNppFile = fullfile('data','processed','npp_carr2002_seawifs_pathfinder_zeuc02.mat');
    fullpathFigureName = 'fig_monthly_npp_carr2002_seawifs_zeuc02';
elseif ~isSensorAquaModis && ~isZeuCarr
    fullpathOutputNppFile = fullfile('data','processed','npp_carr2002_seawifs_pathfinder_zeub97.mat');
    fullpathFigureName = 'fig_monthly_npp_carr2002_seawifs_zeub97';
end

% Query points for interpolation
qLats = (-89.5:1:89.5)';
qLons = (-179.5:1:179.5)';

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD INPUT DATA AND REGRID TO COMMON GRID
% -------------------------------------------------------------------------

% Load input data
load(fullpathInputChla,'chla','chla_lat','chla_lon');
load(fullpathInputPar0,'par0','par0_lat','par0_lon') ;
load(fullpathInputSst,'sst','sst_lat','sst_lon');

% Data grid
[Xc, Yc, Tc] = ndgrid(chla_lat, chla_lon, (1:12)');
[Xp, Yp, Tp] = ndgrid(par0_lat, par0_lon, (1:12)');
[Xs, Ys, Ts] = ndgrid(sst_lat, sst_lon, (1:12)');

% Query grid
[qX, qY, qT] = ndgrid(qLats, qLons, (1:12)');

% Interpolant -use first-order (linear) interpolation and extrapolation
Fchla = griddedInterpolant(Xc, Yc, Tc, chla, 'linear'); 
Fpar0 = griddedInterpolant(Xp, Yp, Tp, par0, 'linear');
Fsst  = griddedInterpolant(Xs, Ys, Ts, sst, 'linear'); 

% Regrid chla, PAR0 and SST from the default grids used by AVHRR 
% Pathfinder and SeaWiFS to the grid defined by qLats and qLons. 
% Regridding is the process of interpolating from one grid resolution 
% to a different one. 
qChla = Fchla(qX, qY, qT);
qPar0 = Fpar0(qX, qY, qT);
qSst = Fsst(qX, qY, qT);
% figure(); pcolor(flipud(rot90(qSst(:,:,1)))); caxis([-2 25]); shading interp;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - APPLY CARR 2002 ALGORITHM
% -------------------------------------------------------------------------

qNppCarr = Carr2002algorithm(qChla,qSst,qPar0,isZeuCarr); % mg C m-2 d-1

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - INSPECTION AND SAVE
% -------------------------------------------------------------------------

% Save output
npp_lat = qLats;
npp_lon = qLons;
npp_avg = qNppCarr;
save(fullpathOutputNppFile,'npp_avg','npp_lat','npp_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputNppFile,[],'mg C m^{-2} d^{-1}',...
    0,1000,true,fullpathFigureName,'NPP Carr 2002')
