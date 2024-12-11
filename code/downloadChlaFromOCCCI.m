function downloadChlaFromOCCCI(fullpathInputChlaOccciDir)

% DOWNLOADCHLAFROMOCCI This function programmatically downloads chla data 
% from OC-CCI.
%
%   INPUT: 
%       fullpathInputChlaOccciDir - directory to store the downloaded .nc files
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

%% Initialise

% Initialise the common parts of the URL
commonSuffix = "https://www.oceancolour.org/thredds/ncss/cci/v6.0-release/geographic/monthly/chlor_a/";
commonMiddle = "?var=chlor_a&disableLLSubset=on&disableProjSubset=on&horizStride=1&";

% Initialise the start year and end year
startYear = 1997;
endYear = 2024;
endMonth = 6; % stop at June for 2024

%% Loop over the years

for iYear = startYear:endYear

    % Determine the start and end months for each year
    if iYear == startYear
        startMonth = 9; % start from September 1997
    else
        startMonth = 1; % for other years, start from January
    end
    if iYear == endYear
        lastMonth = endMonth; % stop at June 2024
    else
        lastMonth = 12; % for other years, go until December
    end

    % Loop over the months for the current year
    for iMonth = startMonth:lastMonth
        monthStr = sprintf('%02d', iMonth);

        if (iYear == startYear && iMonth == startMonth)
            dateStr = sprintf("%04d-%02d-04T00%%3A00%%3A00Z", iYear, iMonth);
        else
            dateStr = sprintf("%04d-%02d-01T00%%3A00%%3A00Z", iYear, iMonth);
        end

        generatedURL = strcat(commonSuffix, num2str(iYear), ...
            "/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-", ...
            num2str(iYear), monthStr, "-fv6.0.nc", commonMiddle, ...
            "time_start=", dateStr, "&time_end=", dateStr, "&timeStride=1&accept=netcdf");

        % Output file name 
        filename = strcat("/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-",...
            num2str(iYear), monthStr, "-fv6.0.nc");

        % Download the file
        websave(strcat(fullpathInputChlaOccciDir,filename), generatedURL);

    end % iMonth
end % iYear

end