
% ======================================================================= %
%                                                                         %
% This script visually compares two interpolation methods (standard and   %
% custom) on a global chlorophyll a concentration dataset (OC-CCI).       %
%                                                                         % 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
%   Version 1.1 - 11 Dec 2024 (improved code)                             %
%                                                                         %
% ======================================================================= %

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

% Output files
fullpathOutputInterpFile1 = fullfile('data','processed','interp_custom_chla_occci.mat');
fullpathOutputInterpFile2 = fullfile('data','processed','interp_standard_chla_occci.mat');

% Load the dataset
load(fullfile('data','processed','chla_occci.mat'),'chla','chla_lat','chla_lon')
lats = chla_lat;
lons = chla_lon;
nLats = length(chla_lat);
nLons = length(chla_lon);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - CUSTOM INTERPOLATION METHOD
% -------------------------------------------------------------------------

chlaGridOriginal = chla;

% Interpolate in time using a 1D time slice vector, improved with respect
% to default interp1
chlaInterpolatedMethod1 = timeInterpolationCustom(chlaGridOriginal,(1:12));

% Save
chla = chlaInterpolatedMethod1;
save(fullpathOutputInterpFile1,'chla','chla_lat','chla_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputInterpFile1,[],'mg m^{-3}',...
    0,1,true,'fig_interp_custom_chla_occci','Chla OC-CCI interp custom')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - STANDARD INTERPOLATION METHOD
% -------------------------------------------------------------------------

chlaGridOriginal = chla;

% Interpolate in time using a 1D time slice vector, standard interp1
chlaInterpolatedMethod2 = timeInterpolationStandard(chlaGridOriginal,(1:12)); % standard time interpolation method

% Save
chla = chlaInterpolatedMethod2;
save(fullpathOutputInterpFile2,'chla','chla_lat','chla_lon','-v7.3')

% Visual inspection
prepareDataForPlotting(fullpathOutputInterpFile2,[],'mg m^{-3}',...
    0,1,true,'fig_interp_standard_chla_occci','Chla OC-CCI interp standard')
