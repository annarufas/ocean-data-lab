function qArray = standardTimeInterpolation(array,t)

% STANDARDTIMEINTERPOLATION Fills temporal gaps in a 3D array using 
% interpolation. This function fills temporal gaps (NaN values) in a 3D 
% array of dimensions lat x lon x time by interpolating over the time 
% dimension. Interpolation is applied only for grid cells with at least 2 
% non-zero values (i.e., 2 months of values). NaN values in areas without 
% sufficient data or invalid land regions are retained.
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
             missingDataMask = (timeSeries==0);
             qArray(iRow,iCol,missingDataMask) = interp1(t(~missingDataMask),...
                 timeSeries(~missingDataMask),t(missingDataMask),'nearest','extrap')';
         end
    end
end

%% Post-processing: set land pixels or those with less than 2 months of data to NaN 

qArray(qArray == 0) = NaN;

end % standardTimeInterpolation