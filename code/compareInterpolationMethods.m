
% ======================================================================= %
%                                                                         %
% This script visually compares three interpolation methods on a global   %
% chlorophyll a concentration dataset.                                    %
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

% Output files
fullpathOutputInterpFile1 = fullfile('figures','interp_m1_chla_occci.mat');
fullpathOutputInterpFile2 = fullfile('figures','interp_m2_chla_occci.mat');
fullpathOutputInterpFile3 = fullfile('figures','interp_m3_chla_occci.mat');

% Load the dataset
load(fullfile('data','processed','chla_occci.mat'),'chla','chla_lat','chla_lon')
lats = chla_lat;
lons = chla_lon;
nLats = length(chla_lat);
nLons = length(chla_lon);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - INTERPOLATION METHOD 1
% -------------------------------------------------------------------------

original_var = chla;

% Set to 0 pixels with NaN as we need to interpolate in such pixels
original_var(isnan(original_var)) = 0;

% Interpolate in time using a 1D time slice vector, improved with respect
% to default interp1
[interp_var] = cleverTimeInterpolation(original_var,[1:12]',lats,lons);

% Save
chla = interp_var;
save(fullpathOutputInterpFile1,'chla','chla_lat','chla_lon','-v7.3')

% Visual inspection
plotMonthlyMaps(fullpathOutputInterpFile1,[],'mg m^{-3}',...
    0,1,true,[],'fig_interp_m1_chla_occci','Chla OC-CCI interp M1')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - INTERPOLATION METHOD 2
% -------------------------------------------------------------------------

original_var = chla;

% Set to 0 pixels with NaN as we need to interpolate in such pixels
original_var(isnan(original_var)) = 0;

% Interpolate in time using a 1D time slice vector, standard interp1
interp_var = chla;
for iLon = 1:nLons   
    for iLat = 1:nLats
         slicetime = squeeze(original_var(iLon,iLat,:));
         % If in this polar region there have been >= 2 months when there 
         % was data, interpolate and extrapolate (we need at least 2 points 
         % to extrapolate)
         if (sum(slicetime > 0) > 1) 
             [t] = (1:12);
             matchZeros = (slicetime==0);
             interp_var(iLon,iLat,matchZeros) = interp1(t(~matchZeros),...
                 slicetime(~matchZeros),t(matchZeros),'nearest','extrap')';
         end
    end
end

% Set to NaN land areas / areas that have been covered by clouds for >= 11 months
interp_var(interp_var==0) = NaN; 

% Save
chla = interp_var;
save(fullpathOutputInterpFile2,'chla','chla_lat','chla_lon','-v7.3')

% Visual inspection
plotMonthlyMaps(fullpathOutputInterpFile2,[],'mg m^{-3}',...
    0,1,true,[],'fig_interp_m2_chla_occci','Chla OC-CCI interp M2')

