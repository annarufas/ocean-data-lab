
% ======================================================================= %
%                                                                         %
%               Zeu climatology from CMEMS's kd product                   % 
%                                                                         %
% This script creates a global gridded climatology of euphotic layer      %
% depth (zeu) using the attenuation light coefficient (kd) product from   %
% Copernicus Marine Service (CMEMS). The output array is 1080 x 2160      % 
% pixels and has units of m.                                              %                                        
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
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

% Output directory
fullpathOutputZeuFile = fullfile('data','processed','zeu_calculated_onepercentpar0.mat');

% Load kd data
fullpathKdClimatology = fullfile('data','processed','kd_cmems_bgc.mat');
load(fullpathKdClimatology,'kd','kd_lat','kd_lon') % m-1

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - CALCULATIONS
% -------------------------------------------------------------------------

% Standard formula to calculate zeu as the depth where PAR is 1% of PAR0,
% as defined in the publication of Buesseler et al. (2020)
% (www.pnas.org/cgi/doi/10.1073/pnas.1918114117)

zeu = -log(0.01)./kd;     

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - INSPECTION AND SAVE
% -------------------------------------------------------------------------

% Save output
zeu_lat = kd_lat;
zeu_lon = kd_lon;
save(fullfile(fullpathOutputZeuFile),'zeu','zeu_lat','zeu_lon')   

% Visual inspection
plotMonthlyMaps(fullpathOutputZeuFile,[],'m',...
    0,200,true,[],'fig_monthly_zeu_calculated_onepercentpar0','Zeu calculated from kd')
