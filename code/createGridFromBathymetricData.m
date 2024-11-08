
% ======================================================================= %
%                                                                         %
%                        Grid from GEBCO bathymetry                       % 
%                                                                         %
% This script creates a grid from the General Bathymetric Chart of the    %
% Ocean (GEBCO) bathymetry. This grid is used in subsequent analyses,     %
% such as sampling ocean locations or regridding other oceanographic      %
% datasets. A grid of high resolution and another of lower resolution are % 
% created.                                                                %        
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
%                                                                         %
% ======================================================================= %

% Clear workspace, close figures
close all; clear all; clc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

fullpathInputBathymetricData = fullfile('data','processed','bathymetry_gebco.mat');
fullpathOutputGridFileHighRes = fullfile('data','raw','grid','grid_GEBCO_2160_1080.mat');
fullpathOutputGridFileLowRes = fullfile('data','raw','grid','grid_GEBCO_360_180.mat');

% Parameters
maxSeafloorDepth = 5000; % set maximum depth threshold at 5000 meters (optimises storage)
depthInterval = 10;      % depth interval for grid layers

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT BATHYMETRIC DATA
% -------------------------------------------------------------------------

% Load bathymetry
load(fullpathInputBathymetricData,'bathym','bathym_lat','bathym_lon'); 
latHighRes = bathym_lat;
lonHighRes = bathym_lon;

% Matrix for storing seafloor depths (high resolution)
Dhigh = -1.*bathym; % transform into non-negative numbers
Dhigh(isnan(Dhigh)) = 0;
Dhigh(Dhigh > maxSeafloorDepth) = maxSeafloorDepth; % defie maximum depth limit

% Create lower resolution matrix
latLowRes = (-89.5:1:89.5);
lonLowRes = (-179.5:1:179.5);
[X,Y] = ndgrid(bathym_lat,bathym_lon); % original array
[qX,qY] = ndgrid(latLowRes,lonLowRes); % query points for interpolation 
F = griddedInterpolant(X,Y,Dhigh);
Dlow = F(qX,qY); 

figure()
pcolor(Dlow)
shading interp
colormap(jet)
box on
colorbar
title('Seafloor depths')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - GENERATE THE GRID AND SAVE IT
% -------------------------------------------------------------------------

[x,y,ixBb,iyBb,Xbb,Ybb,Zsb3d,Zsb] = generateGrid(Dhigh,maxSeafloorDepth,...
    depthInterval,lonHighRes,latHighRes);
save(fullpathOutputGridFileHighRes,'x','y','ixBb','iyBb','Xbb','Ybb','Zsb3d','Zsb','-v7.3');
clear x y ixBb iyBb Xbb Ybb Zsb3d Zsb

[x,y,ixBb,iyBb,Xbb,Ybb,Zsb3d,Zsb] = generateGrid(Dlow,maxSeafloorDepth,...
    depthInterval,lonLowRes,latLowRes);
save(fullpathOutputGridFileLowRes,'x','y','ixBb','iyBb','Xbb','Ybb','Zsb3d','Zsb','-v7.3');
clear x y ixBb iyBb Xbb Ybb Zsb3d Zsb
 
% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTION USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function [x,y,ixBb,iyBb,Xbb,Ybb,Zsb3d,Zsb] = generateGrid(D,maxSeafloorDepth,...
    depthInterval,lons,lats)

    nx = length(lats);
    ny = length(lons);
    x = lats;
    y = lons;

    nz = maxSeafloorDepth/depthInterval; % number of depth layers
    Zsb3d = NaN(nx,ny,nz);               % 3D matrix for depth layers
    Zsb = NaN(nz,nx*ny);                 % 2D matrix for depth layers

    % Initialise arrays for storing non-empty grid cells
    ixBb = zeros(nx*ny, 1);
    iyBb = zeros(nx*ny, 1);
    Xbb = zeros(nx*ny, 1);
    Ybb = zeros(nx*ny, 1);

    iLoc = 0; % counter for non-empty cells
    
    for iLat = 1:nx
        for iLon = 1:ny
            localDepths = (depthInterval/2 : depthInterval : D(iLat,iLon));
            if ~isempty(localDepths)

                % Update 3D and 2D matrices with depth values
                Zsb(1:length(localDepths),iLoc+1) = localDepths;
                Zsb3d(iLat,iLon,1:length(localDepths)) = localDepths;

                % Store coordinates for non-empty cells
                iLoc = iLoc + 1;
                ixBb(iLoc) = iLat;
                iyBb(iLoc) = iLon;
                Xbb(iLoc) = x(iLat);
                Ybb(iLoc) = y(iLon);
                
            end  
        end
    end

    % Trim empty elements from coordinate arrays
    ixBb = ixBb(1:iLoc);
    iyBb = iyBb(1:iLoc);
    Xbb = Xbb(1:iLoc);
    Ybb = Ybb(1:iLoc);
    
end % generateGrid
