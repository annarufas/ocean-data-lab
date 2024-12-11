function qArray = cleverTimeInterpolation(array,t)

% CLEVERTIMEINTERPOLATION Fills temporal gaps in a 3D array using multiple 
% interpolation strategies. This function performs multi-step interpolation 
% to fill gaps in a 3D array (lat x lon x time). It applies three 
% interpolation rounds: (1) interpolation for isolated gaps (zeros between 
% valid data), (2) interpolation for chunks multiple consecutive zeros 
% bounded by valid data, and (3) interpolation for cases where zeros are at 
% the beginning or end of the series, with valid data only on one side 
% (uses flip-based interpolation to extrapolate values from the opposite end).
%
%   INPUT: 
%       array  - 3D array (lat x lon x time) containing gaps (NaN values)
%       t      - 1D time vector (indices) corresponding to the 3rd dimension of array
%
%   OUTPUT:
%       qArray - 3D array (lat x lon x time) with gaps filled using interpolation
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 9 Dec 2024  
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Input validation

if ndims(array) ~= 3
    error('Input array must be a 3D matrix of dimensions lat x lon x time.');
end
if length(t) ~= size(array, 3)
    error('Time vector length must match the 3rd dimension of the array.');
end

%% Replace NaN with 0 before interpolation

array(isnan(array)) = 0; 

%% Interpolation

qArray = array;

for iRow = 1:size(array,1)
    for iCol = 1:size(array,2)
        
        % Extract the time series for the current grid cell
        timeSeries = squeeze(array(iRow,iCol,:));  

        % Perform interpolation only if at least two non-zero values exist
        if sum(timeSeries > 0) > 1 
             %%
            % First round: interpolate isolated gaps
            timeSeries = interpolateIsolatedGaps(timeSeries, t);
%%
            % Second round: interpolate consecutive gaps bounded by valid data
            timeSeries = interpolateBoundedGaps(timeSeries, t);
%%
            % Third round: handle gaps at the edges using flip-based ("capicúa") interpolation
            timeSeries = interpolateEdgeGaps(timeSeries);
 
            % Store the interpolated time series
            qArray(iRow,iCol,:) = timeSeries;

         else 

            % If insufficient data, retain original values
            qArray(iRow,iCol,:) = timeSeries;

         end
    end
end

%% Post-processing: set negative values to zero and zeros to NaN (land)

qArray(qArray<0) = 0;
qArray(qArray==0) = NaN; 

end % cleverTimeInterpolation

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function timeSeries = interpolateIsolatedGaps(timeSeries, t)
% Handles isolated gaps (0 values) between valid data points.
    
    missingDataMask = (timeSeries == 0);
    if any(missingDataMask)
        % Interpolate only for missing values surrounded by valid data
        for i = 2:length(t) - 1
            if timeSeries(i) == 0 && timeSeries(i - 1) > 0 && timeSeries(i + 1) > 0
                timeSeries(i) = interp1(t(~missingDataMask),...
                    timeSeries(~missingDataMask),t(i),'linear');
            end
        end
    end
    
end % interpolateIsolatedGaps

% *************************************************************************

function timeSeries = interpolateBoundedGaps(timeSeries, t)
% Handles consecutive gaps bounded by valid data on both sides.

    for idxAttempt = 1:4 % repeat up to 4 times for multiple chunks
        missingDataMask = (timeSeries == 0);
        if ~any(missingDataMask)
            break; % exit if no gaps remain
        end

        % Locate the start and end of each consecutive gap
        gapStart = find(diff([0; missingDataMask; 0]) == 1);
        gapEnd = find(diff([0; missingDataMask; 0]) == -1) - 1;

        for i = 1:length(gapStart)
            % Interpolate values within the gap
            if gapStart(i) > 1 && gapEnd(i) < length(t)
                validPoints = [gapStart(i) - 1, gapEnd(i) + 1];
                timePoints = [t(validPoints(1)), t(validPoints(2))];
                interpValues = [timeSeries(validPoints(1)), timeSeries(validPoints(2))];
                timeSeries(gapStart(i):gapEnd(i)) = interp1(timePoints, interpValues, t(gapStart(i):gapEnd(i)), 'linear');
            end
        end
    end
    
end % interpolateBoundedGaps

% *************************************************************************

function timeSeries = interpolateEdgeGaps(timeSeries)
% Handles gaps at the edges of the time series using flip-based extrapolation.

    missingDataMask = (timeSeries == 0);
    if any(missingDataMask)
        
        nGaps = sum(missingDataMask);
        
        % Identify the first and last valid values
        firstValid = find(timeSeries > 0,1,'first');
        lastValid = find(timeSeries > 0,1,'last');
        
        % Prepare interpolation vector, identify zero entries in the interpolation vector
        wrapGapTimeSeries = zeros(nGaps+2, 1);
        wrapGapTimeSeries(1) = timeSeries(firstValid);
        wrapGapTimeSeries(end) = timeSeries(lastValid);
        matchZerosVec = (wrapGapTimeSeries==0);
        wrapGapTime = 1:nGaps+2;

        % Vector filled in
        gapFilled = interp1(wrapGapTime(~matchZerosVec),wrapGapTimeSeries(~matchZerosVec),...
            wrapGapTime(matchZerosVec),'linear','extrap');

        % Three cases
        if (firstValid > 1 && lastValid < 12) % ------#####-----
            
            % Fill both ends of timeSeries
            nValsFirstSection = firstValid-1;
            timeSeries(1:nValsFirstSection) = flip(gapFilled(1:nValsFirstSection));
            timeSeries(lastValid+1:12) = flip(gapFilled(nValsFirstSection+1:end));
         
        elseif (firstValid == 1 && lastValid < 12) % #######--------
        
            % Fill the end of timeSeries
            timeSeries(lastValid+1:12) = flip(gapFilled(:)); 

        elseif (firstValid > 1 && lastValid == 12) % -------#######
            
            % Fill the start of timeSeries
            timeSeries(1:firstValid-1) = flip(gapFilled(:)); 

        end
    end       
             
end % interpolateEdgeGaps

% *************************************************************************
