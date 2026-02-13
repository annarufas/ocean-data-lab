function [zeu,kpar] = calculateZeuFromChlaAndMLD(chla,mld,isOnePercentPar0)

    % Equations as in Fox et al. (2024), after Morel & Maritorena (2001)
    % and Morel et al. (2007)

    % Preallocate arrays
    kd   = NaN(size(chla));
    kpar = NaN(size(chla));
    zeu  = NaN(size(chla));

    % Identify valid input indices
    validIdx = ~isnan(chla) & ~isnan(mld); % & chla > 0;

    % Calculate diffuse attenuation coefficient at 490 nm
    kd(validIdx) = 0.0166 + 0.07242.*(chla(validIdx).^0.68955);

    % Secondary valid mask: only calculate kpar where kd is also valid
    validKd = validIdx & kd > eps;  % eps avoids tiny numbers ≈ 0

    % Compute 1./kd only where valid
    invKd = NaN(size(kd));
    invKd(validKd) = 1 ./ kd(validKd);

    % Define shallow and deep MLD relative to 1/kd
    shallowIdx = validKd & mld <= invKd;
    deepIdx    = validKd & mld >  invKd;

    % Calculate kpar with respective formulas
    kpar(shallowIdx) = 0.0864 + 0.884 * kd(shallowIdx) - (0.00137 ./ kd(shallowIdx));
    kpar(deepIdx)    = 0.0665 + 0.874 * kd(deepIdx)    - (0.00121 ./ kd(deepIdx));

    % Calculate zeu only where kpar is valid
    % Define zeu based on light penetration (1% or 0.1%), as defined in the 
    % publication of Buesseler et al. (2020)
    % (https://www.pnas.org/doi/full/10.1073/pnas.1918114117)
    validKpar = ~isnan(kpar) & kpar > 0;
    if isOnePercentPar0
        zeu(validKpar) = -log(0.01) ./ kpar(validKpar);  % 1% light level (= 4.6 / kpar)
    else
        zeu(validKpar) = -log(0.001) ./ kpar(validKpar); % 0.1% light level (= 6.9 / kpar)
    end

end % calculateZeuFromChlaAndMLD