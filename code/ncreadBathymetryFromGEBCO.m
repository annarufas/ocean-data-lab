
% ======================================================================= %
%                                                                         %
%                             GEBCO bathymetry                            %
%                                                                         %
% This script processes bathymetric data from the General Bathymetric     %
% Chart of the Ocean (GEBCO). I have chosen the file GEBCO_2024 Grid      %
% (sub-ice topo/bathy). The high-resolution dataset (43200 x 86400) is    %
% downsampled to a manageable 1080 x 2160 grid, reducing resolution x40   %
% for optimised computational efficiency in further analyses.             %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global
%   - Space resolution: 86400 x 43200                                     %                     
%   - Version: 2024                                                       %
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

filenameInputGebco = 'GEBCO_2024_sub_ice_topo.nc';
filenameOutputGebco = 'bathymetry_gebco.mat';
fullpathInputGebcoFile = fullfile('data','raw','bathymetry',filenameInputGebco);
fullpathOutputGebcoFile = fullfile('data','processed',filenameOutputGebco);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT THE DATA
% -------------------------------------------------------------------------

S = ncinfo(fullpathInputGebcoFile); % short summary

% Read in data
lon = single(ncread(fullpathInputGebcoFile,'lon'));
lat = single(ncread(fullpathInputGebcoFile,'lat'));
bathymetry = single(ncread(fullpathInputGebcoFile,'elevation'));

% Some arrangements
bathymetry(bathymetry > 0) = NaN; % land

% Swap lon and lat dimensions to get lat x lon
bathymetry_perm = permute(bathymetry, [2, 1]);

% Regrid to a lower resolution, from 43200 x 86400 to 1080 x 2160 (decrease
% by a factor of 40)
newLat = linspace(min(lat),max(lat),numel(lat)/40)';
newLon = linspace(min(lon),max(lon),numel(lon)/40)';
[X, Y] = ndgrid(lat,lon); % original array
[qX, qY] = ndgrid(newLat,newLon); % query points for interpolation 
F = griddedInterpolant(X, Y, bathymetry_perm);
qBathym = squeeze(F(qX, qY)); 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CHECKINGS AND SAVE
% -------------------------------------------------------------------------

% Visual inspection
figure()
pcolor(qBathym); 
caxis([-7000 0]);
cb = colorbar('FontSize', 12); 
cb.Label.String = 'Depth (m)';
shading interp
colormap(jet)
box on

% Save the data
bathym = qBathym;
bathym_lon = newLon;
bathym_lat = newLat;
save(fullpathOutputGebcoFile,'bathym','bathym_lat','bathym_lon','-v7.3')
