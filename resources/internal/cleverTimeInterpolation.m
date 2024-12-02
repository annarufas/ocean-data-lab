function qArray = cleverTimeInterpolation(array,t)

% CLEVERTIMEINTERPOLATION This function performs a multi-step interpolation 
% on a time vector. It is designed for filling gaps in time-series. It uses 
% three distinct interpolation rounds to handle different types of missing 
% values: (i) cases where zeros are surrounded by valid data, (ii) cases 
% where multiple consecutive zeros are surrounded at both sides by valid 
% data, (iii) cases where zeros are at the beginning or end of the series, 
% with valid data only on one side (uses "capicúa" or flip-based 
% interpolation to extrapolate values from the opposite end).
%
%   INPUT: 
%       array  - array with gaps
%       t      - indices of time vector
%
%   OUTPUT:
%       qArray - array with gaps filled in
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 13 Nov 2024  
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Get rid of NaNs

qArray = array;

for iRow = 1:size(array,1)
    for iCol = 1:size(array,2)
         timeSlice = squeeze(array(iRow,iCol,:));  

         if sum(timeSlice > 0) > 1 % if there have been >= 2 months in time when there was chla, interpolate or extrapolate

             % First interpolation round: interpolate zeros between 
             % non-zero values

             matchZeros = (timeSlice==0);
             if any(matchZeros)
                 for iMonth = 2:11 % start at second month 
                     if (timeSlice(iMonth) == 0 && timeSlice(iMonth-1) ~= 0 && timeSlice(iMonth+1) ~= 0)
                        timeSlice(iMonth) = interp1(t(~matchZeros),timeSlice(~matchZeros),...
                            t(iMonth),'linear','extrap');
                     end
                 end
             end

             % Second interpolation round: interpolate in bounded spaces,
             % where various zero values are bounded in both sides. There
             % can be more than one chunk meeting this condition. Repeat up
             % to four times

             for i = 1:4
                matchZeros = (timeSlice==0);
                if any(matchZeros)

                     iEndBoundedSpace = 12; % default

                     % Search for the beginning and the end of the zero trail
                     for iMonth = 2:11
                         if (timeSlice(iMonth-1) ~= 0 && timeSlice(iMonth) == 0 && timeSlice(iMonth+1) == 0) 
                             iStartBoundedSpace = iMonth;
                             iEndBoundedSpace = iStartBoundedSpace;
                             for iBound = iStartBoundedSpace+1:12
                                 if (timeSlice(iBound) == 0)
                                     iEndBoundedSpace = iEndBoundedSpace + 1;
                                 else
                                     break
                                 end
                             end
                             break
                         end 
                     end

                     % If the end is not in the last position
                     if (iEndBoundedSpace < 12)
                         nValsToInterp = iEndBoundedSpace - iStartBoundedSpace + 1;
                         limStart = timeSlice(iStartBoundedSpace-1);
                         limEnd = timeSlice(iEndBoundedSpace+1);
                         vectorToInterp = zeros(nValsToInterp+2,1);
                         vectorToInterp(1) = limStart;
                         vectorToInterp(end) = limEnd;
                         matchZerosVec = (vectorToInterp==0);
                         tt = 1:nValsToInterp+2;
                         interpVals = interp1(tt(~matchZerosVec),vectorToInterp(~matchZerosVec),...
                             tt(matchZerosVec),'linear','extrap');
                         timeSlice(iStartBoundedSpace:iEndBoundedSpace) = interpVals;
                     end

                else

                    break

                end

             end % number of chunks for repetition

             % Third interpolation round: capicua    
             matchZeros = (timeSlice==0);
             if any(matchZeros)

                 nValsToInterp = sum(matchZeros);
                 iStartMatchZero = find(timeSlice == 0, 1,'first');
                 iEndMatchZero = find(timeSlice == 0,1,'last');

                 if (iStartMatchZero > 1 && iEndMatchZero < 12) % #####------##### (# = val)
                     limStart = timeSlice(iStartMatchZero-1);
                     limEnd = timeSlice(iEndMatchZero+1);
                     vectorToInterp = zeros(nValsToInterp+2,1);
                     vectorToInterp(1) = limStart;
                     vectorToInterp(end) = limEnd;
                     matchZerosVec = (vectorToInterp==0);
                     tt = 1:nValsToInterp+2;
                     interpVals = interp1(tt(~matchZerosVec),vectorToInterp(~matchZerosVec),...
                         tt(matchZerosVec),'linear','extrap');
                     timeSlice(iStartMatchZero:iEndMatchZero) = interpVals;
                     qArray(iRow,iCol,:) = timeSlice; 
                 else
                     iLimStart = find(timeSlice(:),1,'first');
                     iLimEnd = find(timeSlice(:),1,'last');
                     limStart = timeSlice(iLimStart);
                     limEnd = timeSlice(iLimEnd);
                     vectorToInterp = zeros(nValsToInterp+2,1);
                     vectorToInterp(1) = limStart;
                     vectorToInterp(end) = limEnd;
                     matchZerosVec = (vectorToInterp==0);
                     tt = 1:nValsToInterp+2;
                     interpVals = interp1(tt(~matchZerosVec),vectorToInterp(~matchZerosVec),...
                         tt(matchZerosVec),'linear','extrap');
                     % Three cases
                     if (iLimStart > 1 && iLimEnd < 12) % ------#####-----
                         nValsFirstSec = iLimStart-1;
                         timeSlice(1:iLimStart-1) = flip(interpVals(1:nValsFirstSec));
                         timeSlice(iLimEnd+1:12) = flip(interpVals(nValsFirstSec+1:nValsToInterp));
                         qArray(iRow,iCol,:) = timeSlice; 
                     elseif (iLimStart == 1 && iLimEnd < 12) % #######--------
                         timeSlice(iLimEnd+1:12) = flip(interpVals(:)); 
                         qArray(iRow,iCol,:) = timeSlice; 
                     elseif (iLimStart > 1 && iLimEnd == 12) % -------#######
                         timeSlice(1:iLimStart-1) = flip(interpVals(:)); 
                         qArray(iRow,iCol,:) = timeSlice; 
                     end
                 end

             else 

                 qArray(iRow,iCol,:) = timeSlice;

             end

         else % sum(timeSlice > 0) <= 1

            qArray(iRow,iCol,:) = timeSlice;

         end
    end
end

% Post-processing: set negative values to zero and zeros to NaN (land)
qArray(qArray<0) = 0;
qArray(qArray==0) = NaN; 

end % cleverTimeInterpolation