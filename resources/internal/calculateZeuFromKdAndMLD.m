function [zeu,kpar] = calculateZeuFromKdAndMLD(kd,mld,isOnePercentPar0)

    % Equations as in Fox et al. (2024), after Morel & Maritorena (2001)
    % and Morel et al. (2007)

    % Initialise kpar
    kpar = NaN(size(mld));

    % Calculate diffuse attenuation coefficient for Photosynthetically Available Radiation for MLD <= 1/kd
    kpar(mld <= (1 ./ kd)) = 0.0864 + 0.884 * kd(mld <= (1 ./ kd)) - (0.00137 ./ kd(mld <= (1 ./ kd)));
    
    % Calculate diffuse attenuation coefficient for Photosynthetically Available Radiation for MLD > 1/kd
    kpar(mld > (1 ./ kd)) = 0.0665 + 0.874 * kd(mld > (1 ./ kd)) - (0.00121 ./ kd(mld > (1 ./ kd)));
    
    % Define zeu based on light penetration (1% or 0.1%), as defined in the 
    % publication of Buesseler et al. (2020)
    % (https://www.pnas.org/doi/full/10.1073/pnas.1918114117)
    
    if isOnePercentPar0
        zeu = -log(0.01) ./ kpar;
    else
        zeu = -log(0.001) ./ kpar; 
    end
    
end % calculateZeuFromKdAndMLD
