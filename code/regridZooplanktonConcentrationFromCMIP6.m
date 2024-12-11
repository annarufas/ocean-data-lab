function regridZooplanktonConcentrationFromCMIP6(filenameInput,filenameOutput,...
    filenameGridDomain)

% REGRIDZOOPLANKTONCONCENTRATIONFROMCMIP6 Regrids from native to regular
% grid datasets of zooplankton concentration. This function is run in a
% shell script in SLURM.
%
%   INPUT: 
%       filenameInput      - name of the file containing the zooplankton array in native grid
%       filenameOutput     - name of the file containing the zooplankton array in regular grid
%       filenameGridDomain - grid file to be used for regridding from native to regular
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 29 Oct 2024  
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Load MITgcm ocean model grid

parentDir = fullfile('..', '..'); % go up two directories
fullPathGridDomainDir = fullfile(parentDir,'grid',strcat(filenameGridDomain,'.mat'));
load(fullPathGridDomainDir,'Xbb','Ybb');
        
%% Load arrays to be regridded

load(strcat(filenameInput,'.mat'),'mesozooClimatologyNative','lats','lons')

%% Regrid

% Query points for interpolation
lonRegular = Xbb(:);
latRegular = Ybb(:);

% If the climatology is 4D, it has depth dimensions and thus the data are in 
% units of concentration (mg C m-3)
if (ndims(mesozooClimatologyNative) == 4)
    nDepths = size(mesozooClimatologyNative,3);
    mesozooClimatologyRegular = NaN(length(Xbb),nDepths,12); 
    for iDepth = 1:nDepths
        mesozooClimatologyRegular(:,iDepth,:) = interp_2dfield(squeeze(mesozooClimatologyNative(:,:,iDepth,:)),...
            double(lons),double(lats),lonRegular,latRegular,[],[],0);
    end
% If the climatology is 3D, it has been depth-integrated and thus the data are in 
% units of mg C m-2
else   
    mesozooClimatologyRegular = NaN(length(Xbb),12); 
    mesozooClimatologyRegular(:,:) = interp_2dfield(squeeze(mesozooClimatologyNative(:,:,:)),...
        double(lons),double(lats),lonRegular,latRegular,[],[],0);
end
    
save(strcat(filenameOutput,'.mat'),'mesozooClimatologyRegular');

end