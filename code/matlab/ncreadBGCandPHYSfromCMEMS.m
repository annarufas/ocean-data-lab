
% ======================================================================= %
%                                                                         %
%          Chla, NPP and kd (BGC), MLD, T and ice fraction (PHYS)         %
%                       climatologies from CMEMS                          %
%                                                                         %
% This script reads in Copernicus Marine Service (CMEMS) monthly products %
% of biogeochemical (BGC) and physical (PHYS) variables downloaded using  %
% the Jupyter Notebook "downloadBGCandPHYSfromCMEMS.ipynb". It regrids    %
% the data if downloaded datasets are too large and uses the algorithm    %
% worstcase to compute the climatological values of uncertainty arrays.   %
%                                                                         %
% Dataset characteristics:                                                %
%   - https://data.marine.copernicus.eu/products                          %
%   - Time resolution: monthly (BGC) and monthly climatologies (PHYS)     %
%   - Space resolution: 4 x 4 km (BGC) and 0.083° × 0.083° (PHYS)         %
%   - Version: 2024                                                       %
%   - Units: BGC - chlorophyll a concentration (chla): mg m-3             % 
%            BGC - net primary production (NPP): mg C m-2 d-1             %
%            BGC - diffuse attenuation coefficient (kd): m-1              % 
%            PHYS - mixed layer depth (MLD): m                            %
%            PHYS - seawater temperature: °C                              %
%            PHYS - sea ice area fraction: 0-1 fraction                   %
%                                                                         %                                                                        
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 29 Oct 2024                                   %
%   Version 1.1 - 11 Nov 2024 (added uncertainty array for NPP)           %
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

% As defined in the Jupyter Notebook "downloadBGCandPHYSfromCMEMS.ipynb"
filenameInput = {'mod_bgc_glo_npp.nc',... % 'mod_bgc_glo_npp_uncertainty.nc',...
                 'mod_bgc_glo_chla.nc',...
                 'mod_bgc_glo_kd.nc',...
                 'mod_phys_glo_mld.nc',...
                 'mod_phys_glo_icefrac.nc',...
                 'mod_phys_glo_temp.nc'};

% Outputs contain grouped arrays
filenameOutput = {'npp_cmems_bgc.mat',...
                  'chla_cmems_bgc.mat',...
                  'kd_cmems_bgc.mat',...
                  'mld_cmems_phys.mat',...
                  'icefrac_cmems_phys.mat',...
                  'temp_cmems_phys.mat'};

% Define groups based on the fourth term in filenameOutput
terms = {'npp', 'chla', 'kd', 'mld', 'icefrac', 'temp'};

% Paths to directories
fullpathInputCmemsDir = fullfile('data','raw','CMEMS');
fullpathOutputCmemsDir = fullfile('data','processed');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - DOWNLOAD THE DATA
% -------------------------------------------------------------------------

% Run the Jupyter Notebook "downloadBGCandPHYSfromCMEMS.ipynb" to
% programmatically download data from CMEMS.

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - GROUP THE DATA (MEAN AND UNCERTAINTY)
% -------------------------------------------------------------------------

groupedFiles = struct();

for i = 1:length(filenameInput)
    % Remove the file extension (.nc) and split by '_'
    filenameWithoutExt = strtok(filenameInput{i}, '.');
    parts = strsplit(filenameWithoutExt, '_'); 
    fourthTerm = parts{4}; 
    
    % Check if the group already exists in the struct, if not, initialise it
    if ~isfield(groupedFiles, fourthTerm)
        groupedFiles.(fourthTerm) = {};
    end
    
    % Append the current file to the appropriate group
    groupedFiles.(fourthTerm){end+1} = filenameInput{i};
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - EXTRACT THE DATA
% -------------------------------------------------------------------------

% Naming conventions used by CMEMS
cmemsVarNameLatitude = 'lat';
cmemsVarNameLongitude = 'lon';
cmemsVarNameTime = 'time';
cmemsVarNameDepth = 'depth'; 

% Initialise structure array for output
cmemsData = struct('ID', {}, 'varNames', {}, 'units', {}, 'dataset', {},... 
    'lat', {}, 'lon', {}, 'time', {}, 'depth', {});

% Loop over the grouped files
for iDataset = 1:3 %length(terms)
    term = terms{iDataset};
    files = groupedFiles.(term);
        
    % Initialise the array to store the data for this group
    Dclimall = [];

    % Loop over the files and read them
    for iFile = 1:length(files)
        fileName = files{iFile};
        filePath = fullfile(fullpathInputCmemsDir,fileName);
        
        fprintf('\nReading %s',fileName)
        
        % List variable names and units
        S = ncinfo(filePath); % short summary
        nVars = length(S.Variables);
        varName = cell(nVars,1);
        varUnit = cell(nVars,1);
        for iVar = 1:nVars
            varName{iVar} = S.Variables(iVar).Name;
            iUnit = contains({S.Variables(iVar).Attributes.Name},'units','IgnoreCase',true);
            if (sum(iUnit) > 0) % check that the field 'units' exists before reading it (some variables, like 'flags' don't have units)
                varUnit{iVar} = ncreadatt(filePath,varName{iVar},'units');
            end
        end

        % Get idx of dimensional variable names (different files use different
        % variable names (e.g., some files use 'lat', others use 'latitude'))
        iLat = find(contains(string(varName),cmemsVarNameLatitude,'IgnoreCase',true));
        iLon = find(contains(string(varName),cmemsVarNameLongitude,'IgnoreCase',true));
        iTime = find(contains(string(varName),cmemsVarNameTime,'IgnoreCase',true));
        iDepth = find(contains(string(varName),cmemsVarNameDepth,'IgnoreCase',true));

        % Read longitude, latitude, time and depth (if exists)
        lat = double(ncread(filePath,varName{iLat}));
        lon = double(ncread(filePath,varName{iLon}));
        time = double(ncread(filePath,varName{iTime})); 
        if ~isempty(iDepth)
            depth = double(ncread(filePath,varName{iDepth})); 
        end

        % Transform time into a more conventional format
        timeCalendar = convertTimeUnits(varUnit,time,iTime);

        % Get idx of non-structural variables (i.e., ocean colour information)
        isVarOceanColour = ~ismember(1:nVars, [iLat, iLon, iTime, iDepth]);
        oceanColourVarName = varName(isVarOceanColour);

        % Create 'start' and 'count' arguments for the ncread function. For
        % that, I need to find the arrangement of the dimensional variables in 
        % the ocean colour variables
        for i = 1:nVars
            if (isVarOceanColour(i) == 1)
                iDimLat = find(contains(string({S.Variables(i).Dimensions.Name}),cmemsVarNameLatitude,'IgnoreCase',true)); 
                iDimLon = find(contains(string({S.Variables(i).Dimensions.Name}),cmemsVarNameLongitude,'IgnoreCase',true)); 
                iDimTime = find(contains(string({S.Variables(i).Dimensions.Name}),cmemsVarNameTime,'IgnoreCase',true)); 
                iDimDepth = find(contains(string({S.Variables(i).Dimensions.Name}),cmemsVarNameDepth,'IgnoreCase',true));
                break
            end
        end
        if ~isempty(iDepth)
            startIndices = [1           1           1             1]; 
            newDimCounts = [length(lat) length(lon) length(depth) length(time)]; 
            orderIndices = [iDimLat     iDimLon     iDimDepth     iDimTime];
        else
            startIndices = [1           1           1]; 
            newDimCounts = [length(lat) length(lon) length(time)]; 
            orderIndices = [iDimLat     iDimLon     iDimTime];
        end

        % Sort start and counts according to dimensional arrangement in the
        % original dataset
        [~,sortedIndices] = sort(orderIndices);
        start = startIndices(sortedIndices); % start indices ordered
        count = newDimCounts(sortedIndices); % dimensional information ordered

        % Read data and save it in a standardised format (following WOA):
        % 1st dimension: latitude
        % 2nd dimension: longitude
        % 3rd dimension: depth
        % 4th dimension: time
        % The ordering follows the 'newDimCounts' array specification.

        % Estimate memory requirements before initialising the output array
        sizeDataset = count;
        totalElements = prod(sizeDataset);
        memoryRequired = totalElements * 4; % in bytes
        memoryRequiredMB = memoryRequired / (1024^2); % Convert to MB

        % If the requested memory does not exceed the maximum array size 
        % preference in MATLAB (16 GB), no need to regrid. To this category
        % belong data that are already a climatology (PHYS variables)
        if (memoryRequiredMB < 16000) 

            Dtmp = double(ncread(filePath,oceanColourVarName{1},start,count));
            if ~isempty(iDepth)
                Dperm = permute(Dtmp,[iDimLat iDimLon iDimDepth iDimTime]);
            elseif isempty(iDepth)
                Dperm = permute(Dtmp,[iDimLat iDimLon iDimTime]);
                depth = 0;
            end
            Dclim = Dperm;

            % If this is an uncertainty file, append it to the data (Dclimall)
            if contains(fileName, 'uncertainty')
                Dclim = (Dclim./100).*Dclimall;
                Dclimall = [Dclimall, Dclim];
            else
                Dclimall = Dclim;
            end

        % If the requested memory exceeds the maximum array size preference in 
        % MATLAB (16 GB), regrid the array. To this category belong data that
        % are not yet a climatology (BGC variables)
        else 

            % Regrid the original array to a lower resolution      
            newLon = linspace(-180, 180, 2160)';
            newLat = linspace(-90, 90, 1080)';

            % Determine lat/lon order before regridding
            lonFirstDim = (iDimLon == 1); 

            % Output array dimensions
            if lonFirstDim
                if isempty(iDepth)
                    D = NaN([numel(newLon),numel(newLat),numel(time)],'double'); 
                else
                    D = NaN([numel(newLon),numel(newLat),numel(depth),numel(time)],'double'); 
                end
            else
                if isempty(iDepth)
                    D = NaN([numel(newLat),numel(newLon),numel(time)],'double'); 
                else
                    D = NaN([numel(newLat),numel(newLon),numel(depth),numel(time)],'double'); 
                end
            end

            % Predefine start and count arrays
            if isempty(iDepth)
                start = ones(1,3); 
                count = [count(1),count(2),1];
            else
                start = ones(1,4);
                count = [count(1),count(2),count(3),1];
            end

            % Regrid 2D (lat/lon) array extracted for each time step
            for iTime = 1:length(time)

                % Update indices for time in start array
                if isempty(iDepth)
                    start(3) = iTime;
                else
                    start(4) = iTime;
                end

                % Extract 2D array for regriding
                Dchunk = squeeze(ncread(filePath,oceanColourVarName{1},start,count));
                if isempty(iDepth)
                    if lonFirstDim
                        [X, Y] = ndgrid(lon, lat); % indexes original grid
                        [qX, qY] = ndgrid(newLon, newLat); % low-resolution grid
                    else
                        [X, Y] = ndgrid(lat, lon); 
                        [qX, qY] = ndgrid(newLat, newLon); 
                    end
                    F = griddedInterpolant(X, Y, Dchunk, 'linear', 'none');
                    D(:,:,iTime) = F(qX, qY);
                else
                    if lonFirstDim
                        [X, Y, Z] = ndgrid(lon, lat, (1:length(depth))'); 
                        [qX, qY, qZ] = ndgrid(newLon, newLat, (1:length(depth))'); 
                    else
                        [X, Y, Z] = ndgrid(lat, lon, (1:length(depth))'); 
                        [qX, qY, qZ] = ndgrid(newLat, newLon, (1:length(depth))'); 
                    end
                    F = griddedInterpolant(X, Y, Z, Dchunk, 'linear', 'none');
                    D(:,:,:,iTime) = F(qX, qY, qZ);
                end

            end % iTime
 
            % Calculate data on a climatological basis
            timeCalendarMonths = month(timeCalendar);

            % Initialise output climatological array
            if isempty(iDepth)
                Dclim = NaN([numel(newLat),numel(newLon),12],'double');
            else
                Dclim = NaN([numel(newLat),numel(newLon),numel(depth),12],'double'); 
            end

            % Create the climatology
            for iMonth = 1:12
                idxMonth = (timeCalendarMonths == iMonth);

                if ~isempty(iDepth)
                    if contains(fileName, 'uncertainty') % calculate propagated error (worst-case)
                        Dtmp = handleUncertainty(Dmeansource,D,idxMonth,'depth');
                    else
                        Dtmp = mean(D(:,:,:,idxMonth),4,'omitnan');
                    end
                    Dperm = permute(Dtmp,[iDimLat iDimLon iDimDepth iDimTime]);
                    Dclim(:,:,:,iMonth) = Dperm; 
                elseif isempty(iDepth)
                    if contains(fileName, 'uncertainty') % calculate propagated error (worst-case)
                        Dtmp = handleUncertainty(Dmeansource,D,idxMonth,'nodepth');
                    else
                        Dtmp = mean(D(:,:,idxMonth),3,'omitnan');
                    end
                    Dperm = permute(Dtmp,[iDimLat iDimLon iDimTime]);
                    Dclim(:,:,iMonth) = Dperm;
                end
                
            end % iMonth
            
            if isempty(iDepth)
                depth = 0;
            end

            clear lat lon
            lon = newLon;
            lat = newLat;

        end % checking memory requirements

        % Check for spurious data points
        %figure(); histogram(Dclim(:),100);
        
        % If this is an uncertainty file, append it to the data (Dclimall)
        if contains(fileName, 'uncertainty')
            Dclim = (Dclim./100).*Dclimall;
            Dclimall = [Dclimall, Dclim]; 
        else
            Dclimall = Dclim;  
        end
        
        Dmeansource = D;

    end % iFile   
    
    % Save output
    saveAndVisualiseOutput(term,Dclimall,lat,lon,depth,...
        fullpathOutputCmemsDir,filenameOutput{iDataset})
    
end % iDataset

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

% *************************************************************************

function timeCalendar = convertTimeUnits(varUnit,time,iTime)

% Function to transform time into a calendar (datetime) format. 

    timeUnits = varUnit{iTime};
    timeUnitsNum = extract(timeUnits,digitsPattern);
    
    % Dictionary mapping of units to MATLAB epoch types. I have checked the 
    % 'epochtype' options in the files I want to download
    epochMap = containers.Map({...
        'days since 0000-01-01', ...
        'days since 1900-01-01', ...
        'days since 1904-01-01', ...
        'seconds since 1970-01-01', ...
        'seconds since 1970-01-01 00:00:00'}, ...
        {'datenum', 'excel', 'excel1904', 'posixtime', 'posixtime'});
    
    % Determine the epoch type and date reference
    if epochMap.isKey(timeUnits)
        epochtype = epochMap(timeUnits);
        timeCalendar = datetime(time(:), 'ConvertFrom', epochtype, ...
            'Epoch', [timeUnitsNum{1} '-' timeUnitsNum{2} '-' timeUnitsNum{3}]);
    elseif strcmp(timeUnits, 'hours since 1950-01-01') % this is not an epochtype in MATLAB, handle hours-based unit separately
        refDate = datetime(1950, 1, 1);
        timeCalendar = refDate + hours(time);
    else
        error('Unsupported time unit: %s', timeUnits);
    end

end % convertTimeUnits

% *************************************************************************

function Dtmp = handleUncertainty(Dmeansource,D,idxMonth,type)

% Function to deal with uncertainty arrays and perform error propagation.

    for iRow = 1:size(Dmeansource,1)
        for iCol = 1:size(Dmeansource,2)
            if strcmp(type,'depth')
                
                for iDepth = 1:size(Dmeansource,3)
                    avgData = squeeze(Dmeansource(iRow, iCol, iDepth, idxMonth));
                    errData = squeeze(D(iRow, iCol, iDepth, idxMonth)); 

                    % Check for NaN values for the current month index
                    idxNonNaN = ~isnan(avgData) & ~isnan(errData);
                    if any(idxNonNaN)
                        Davgloc = squeeze(avgData(idxNonNaN)); 
                        Derrloc = squeeze(errData(idxNonNaN)); 
                        if ~isempty(Derrloc)
                            [~,~,~,f_MID,f_UB] = worstcase(@(Davgloc) mean(Davgloc,'omitnan'), Davgloc, Derrloc);
                            Dtmp(iRow,iCol,iDepth,:) = f_UB - f_MID;
                        end
                    end
                end % iDepth
                
            else
                
                avgData = squeeze(Dmeansource(iRow, iCol, idxMonth));
                errData = squeeze(D(iRow, iCol, idxMonth));

                % Check for NaN values for the current month index
                idxNonNaN = ~isnan(avgData) & ~isnan(errData);
                if any(idxNonNaN)
                    Davgloc = squeeze(avgData(idxNonNaN)); 
                    Derrloc = squeeze(errData(idxNonNaN)); 
                    if ~isempty(Derrloc)
                        [~,~,~,f_MID,f_UB] = worstcase(@(Davgloc) mean(Davgloc,'omitnan'), Davgloc, Derrloc);
                        Dtmp(iRow,iCol,:) = f_UB - f_MID;
                    end
                end
                
            end
        end % iCol
    end % iRow 

end % handleUncertainty
            
% *************************************************************************

function saveAndVisualiseOutput(term,Dclimall,lat,lon,depth,...
    fullpathOutputCmemsDir,filenameOutput)
    
    if strcmp(term,'npp')

        npp_avg = Dclimall(:,:,:,1);
        %npp_err = Dclimall(:,:,:,2);
        npp_lat = lat;
        npp_lon = lon;
        save(fullfile(fullpathOutputCmemsDir,filenameOutput),...
            'npp_avg','npp_lat','npp_lon','-v7.3')

        % Visual inspection
        prepareDataForPlotting(fullfile(fullpathOutputCmemsDir,filenameOutput),[],'mg C m-2 d-1',...
            0,1000,true,'fig_monthly_npp_cmems','NPP CMEMS')

    elseif strcmp(term,'chla')

        chla = Dclimall;
        chla_lat = lat;
        chla_lon = lon;
        save(fullfile(fullpathOutputCmemsDir,filenameOutput),...
            'chla','chla_lat','chla_lon','-v7.3')  

        % Visual inspection
        prepareDataForPlotting(fullfile(fullpathOutputCmemsDir,filenameOutput),[],'mg m^{-3}',...
            0,1,true,'fig_monthly_chla_cmems','Chla CMEMS')

    elseif strcmp(term,'kd')

        kd = Dclimall(:,:,:,1);
        %kd_err = Dclimall(:,:,:,2);
        kd_lat = lat;
        kd_lon = lon;
        save(fullfile(fullpathOutputCmemsDir,filenameOutput),...
            'kd','kd_lat','kd_lon','-v7.3')   

        % Visual inspection
        prepareDataForPlotting(fullfile(fullpathOutputCmemsDir,filenameOutput),[],'m^{-1}',...
            0,0.1,true,'fig_monthly_kd_cmems','kd CMEMS')

    elseif strcmp(term,'mld')

        mld = Dclimall;
        mld_lat = lat;
        mld_lon = lon;
        save(fullfile(fullpathOutputCmemsDir,filenameOutput),...
            'mld','mld_lat','mld_lon','-v7.3')   

        % Visual inspection
        prepareDataForPlotting(fullfile(fullpathOutputCmemsDir,filenameOutput),[],'m',...
            0,200,true,'fig_monthly_mld_cmems','MLD CMEMS')

    elseif strcmp(term,'icefrac')

        icefrac = Dclimall;
        icefrac_lat = lat;
        icefrac_lon = lon;
        save(fullfile(fullpathOutputCmemsDir,filenameOutput),...
            'icefrac','icefrac_lat','icefrac_lon','-v7.3')  

        % Visual inspection
        prepareDataForPlotting(fullfile(fullpathOutputCmemsDir,filenameOutput),[],'unitless',...
            0,0.5,true,'fig_monthly_icefrac_cmems','Ice fraction CMEMS')

    elseif strcmp(term,'temp')

        temp = Dclimall;
        temp_lat = lat;
        temp_lon = lon;
        temp_depth = depth;
        save(fullfile(fullpathOutputCmemsDir,filenameOutput),...
            'temp','temp_lat','temp_lon','temp_depth','-v7.3')  

        % Visual inspection
        prepareDataForPlotting(fullfile(fullpathOutputCmemsDir,filenameOutput),18,'ºC',...
            0,25,true,'fig_monthly_temp_cmems','Temperature at 50 m, CMEMS')

    end
    
end % saveAndVisualiseOutput
    
% *************************************************************************
