function [qarray] = cleverTimeInterpolation(array,t,columnCoords,rowCoords)

qarray = array;

for iCol = 1:length(columnCoords)
%     disp(iLat)
    for iRow = 1:length(rowCoords)  
        
         arraymonth = squeeze(array(iRow,iCol,:));  
         
         if (sum(arraymonth > 0) > 1) % if there have been >= 2 months in time when there was chla, interpolate or extrapolate

             % First interpolation round: interpolate in bounded spaces,
             % where a zero value is bounded in both sides
             
             matchZeros = (arraymonth==0);
             nMatchZeros = sum(matchZeros);
             if (nMatchZeros > 0)
                 for iMonth = 2:11 % start at second month 
                     if (arraymonth(iMonth) == 0 && arraymonth(iMonth-1) ~= 0 && arraymonth(iMonth+1) ~= 0)
                        arraymonth(iMonth) = interp1(t(~matchZeros),arraymonth(~matchZeros),...
                            t(iMonth),'linear','extrap');
                     end
                 end
             end

             % Second interpolation round: interpolate in bounded spaces,
             % where various zero values are bounded in both sides. There
             % can be more than one chunk meeting this condition. Repeat up
             % to four times
 
             for i = 1:4
                 
                matchZeros = (arraymonth==0);
                nMatchZeros = sum(matchZeros);
                
                if (nMatchZeros > 0)

                     iEndBoundedSpace = 12; % default

                     % Search for the beginning and the end of the zero trail
                     for iMonth = 2:11
                         if (arraymonth(iMonth-1) ~= 0 && arraymonth(iMonth) == 0 && arraymonth(iMonth+1) == 0) 
                             iStartBoundedSpace = iMonth;
                             iEndBoundedSpace = iStartBoundedSpace;
                             for iBound = iStartBoundedSpace+1:12
                                 if (arraymonth(iBound) == 0)
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
                         limStart = arraymonth(iStartBoundedSpace-1);
                         limEnd = arraymonth(iEndBoundedSpace+1);
                         vectorToInterp = zeros(nValsToInterp+2,1);
                         vectorToInterp(1) = limStart;
                         vectorToInterp(end) = limEnd;
                         matchZerosVec = (vectorToInterp==0);
                         tt = 1:nValsToInterp+2;
                         interpVals = interp1(tt(~matchZerosVec),vectorToInterp(~matchZerosVec),...
                             tt(matchZerosVec),'linear','extrap');
                         arraymonth(iStartBoundedSpace:iEndBoundedSpace) = interpVals;
                     end

                else
                    
                    break
                    
                end
             
             end % number of chunks for repetition

             % Third interpolation round: capicua    
             matchZeros = (arraymonth==0);
             nMatchZeros = sum(matchZeros);
             
             if (nMatchZeros > 0)

                 nValsToInterp = sum(matchZeros);
                 iStartMatchZero = find(arraymonth == 0, 1,'first');
                 iEndMatchZero = find(arraymonth == 0,1,'last');

                 if (iStartMatchZero > 1 && iEndMatchZero < 12) % #####------##### (# = val)
                     limStart = arraymonth(iStartMatchZero-1);
                     limEnd = arraymonth(iEndMatchZero+1);
                     vectorToInterp = zeros(nValsToInterp+2,1);
                     vectorToInterp(1) = limStart;
                     vectorToInterp(end) = limEnd;
                     matchZerosVec = (vectorToInterp==0);
                     tt = 1:nValsToInterp+2;
                     interpVals = interp1(tt(~matchZerosVec),vectorToInterp(~matchZerosVec),...
                         tt(matchZerosVec),'linear','extrap');
                     arraymonth(iStartMatchZero:iEndMatchZero) = interpVals;
                     qarray(iRow,iCol,:) = arraymonth; 
                 else
                     iLimStart = find(arraymonth(:),1,'first');
                     iLimEnd = find(arraymonth(:),1,'last');
                     limStart = arraymonth(iLimStart);
                     limEnd = arraymonth(iLimEnd);
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
                         arraymonth(1:iLimStart-1) = flip(interpVals(1:nValsFirstSec));
                         arraymonth(iLimEnd+1:12) = flip(interpVals(nValsFirstSec+1:nValsToInterp));
                         qarray(iRow,iCol,:) = arraymonth; 
                     elseif (iLimStart == 1 && iLimEnd < 12) % #######--------
                         arraymonth(iLimEnd+1:12) = flip(interpVals(:)); 
                         qarray(iRow,iCol,:) = arraymonth; 
                     elseif (iLimStart > 1 && iLimEnd == 12) % -------#######
                         arraymonth(1:iLimStart-1) = flip(interpVals(:)); 
                         qarray(iRow,iCol,:) = arraymonth; 
                     end
                 end
                     
             else 
                 
                 qarray(iRow,iCol,:) = arraymonth;
                 
             end
             
         else % sum(arraymonth > 0) <= 1
         
            qarray(iRow,iCol,:) = arraymonth;
         
         end
    end
end

qarray(qarray<0) = 0;
qarray(qarray==0) = NaN; % land

end