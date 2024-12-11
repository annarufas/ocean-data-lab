function [par0photoperiod, lengthPhotoperiod] = calculatePAR0fromTrigonometricEquations(...
    latitude, cloudFrac, iceFrac) 

% CALCULATEPAR0FROMTRIGONOMETRICEQUATIONS Calculates photosynthetic active 
% radiation in the surface ocean (PAR0) for a point location from 
% trigonometric/astronomical equations combined with input data of the 
% fraction of clouds and ice covering the sea.
%
%   INPUT: 
%       latitude  - in degrees north
%       cloudFrac - fraction of clouds in oktas (365 x 1 array)
%       iceFrac   - fraction of ice 0-1 (365 x 1 array)
%
%   OUTPUT:
%       par0photoperiod   - in W m-2 (= J s-1 m-2), average PAR0 received during 
%                           the daylight period of each day, 365 x 1 array
%       lengthPhotoperiod - duration (number of hours) of the daylight
%                           period, 365 x 1 array
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

%% Parameter declarations

% Physical constants used in this function
SOLAR_CNT = 1368; % W m-2, average intensity per square meter at the surface of the Earth
EARTH_TILT_FROM_ORBIT_PLANE = 23.45; % deg
ATMOSPHERIC_TRANSMISSION_COEFF = 0.7;
ABSORPTION_BY_VAPOUR_AND_OZONE = 0.09;
PAR_FRAC = 0.43; % fraction of photosynthetically active radiation (PAR)
SURFACE_ALBEDO = 0.06; % standard value for the ocean

% Conversion factors
DEG_TO_RAD = 2*pi/360;
RAD_TO_DEG = 360/(2*pi);
WATTPERM2_TO_UMOLPERM2PERSEC = 1e6/(3.6164e-19*6.02e23); % from W m-2 to umol light m-2 s-1

% [W] = [J s-1]
% E = h x c / lambda 
%   E = energy of a photon (J), h = Plank constant (J s), c = speed of light (m s-1), lambda = photon wavelength (m)
% PAR is between 400-700 nm
% For a photon at 550 nm
% E550 = 6.63e-34 J·s x 3e8 m·s-1 / (500e-9 m) = 3.6164e-19 J per photon
% Avogadro cnt: 1 mol of light = 6.02e23 photons
% 1 mol photons = 1 mol quanta

%% Calculations start here

% Convert
latitudeRad = latitude*DEG_TO_RAD;

% Initialise output arrays
par0photoperiod = zeros(365,1);
lengthPhotoperiod = zeros(365,1);

for iJulianDay = 1:365

    % .....................................................................

    % Irradiance at the top of the atmosphere (Brock 1981)

    % Sun declination angle, the angle between the equatorial plane of the
    % Earth and the axis joining the centres of the Earth and the Sun
    sunDeclinAngleDeg = EARTH_TILT_FROM_ORBIT_PLANE...
        *sin(2*pi*(284+iJulianDay)/365); % between -23.45° and 23.45°
    sunDeclinAngleRad = sunDeclinAngleDeg*DEG_TO_RAD; 
    
    % Zenith angle
    sunZenithAngleRad = acos(sin(sunDeclinAngleRad)*sin(latitudeRad)... 
        + cos(sunDeclinAngleRad)*cos(latitudeRad));

    % Solar noon radiation available at the top of the atmosphere according
    % to latitude
    radiusVector = 1/(sqrt(1+0.033*cos(2*pi*iJulianDay/365))); % between 0.98324 and 1.01671
    irrAtTopAtmosphere = SOLAR_CNT*cos(sunZenithAngleRad)/radiusVector^2; % W m-2

    % .....................................................................

    if (irrAtTopAtmosphere > 0)
        
        % Irradiance at the Earth's surface (Rosati & Miyakoda 1988): only 
        % a percentage of extraterrestrial radiation reaches the surface of
        % the Earth. From this percentage, a fraction reaches the surface 
        % through direct (or beam) radiation and the rest due to diffuse 
        % radiation (atmospheric particles change direction of incident 
        % solar rays)
        
        % Direct component of solar radiation
        irrBeam = irrAtTopAtmosphere...
            *ATMOSPHERIC_TRANSMISSION_COEFF^(1/cos(sunZenithAngleRad)); 
        
        % Diffuse sky radiation under cloudless conditions
        irrDiffuse = ((1-ABSORPTION_BY_VAPOUR_AND_OZONE)*irrAtTopAtmosphere - irrBeam)*0.5; 

        % Total radiation reaching the surface under clear skies at solar
        % noon
        irrAtNoonSurfaceEarthNoCloud = irrBeam + irrDiffuse; % W m-2

        % .....................................................................
        
        % Attenuation by clouds, based on the empirical formula derived by 
        % Reed (1977)

        sunElevationAngleRad = asin(sin(latitudeRad)...
            *sin(EARTH_TILT_FROM_ORBIT_PLANE*DEG_TO_RAD...
            *sin((iJulianDay-82)*DEG_TO_RAD))...
            + cos(latitudeRad)*cos(EARTH_TILT_FROM_ORBIT_PLANE*DEG_TO_RAD...
                *sin((iJulianDay-82)*DEG_TO_RAD)));
        sunElevationAngleDeg = sunElevationAngleRad*RAD_TO_DEG;

        fractionalCloudCover = cloudFrac(iJulianDay)/8; % 8 is the max (oktas)

        irrAtNoonSurfaceEarthCloud = irrAtNoonSurfaceEarthNoCloud ...
            *(1 - 0.62*fractionalCloudCover + 0.0019*sunElevationAngleDeg); % W m-2

        % .................................................................
        
        % Take into account:
        %   albedo effect: a percentage of the radiation that reaches the 
        %       surface of the Earth is reflected (ground reflectance of beam and
        %       diffuse radiation)
        %   ice effect
        %   PAR fraction: only about 43% of the energy of solar radiation is 
        %       actually in the 400 - 700 nm (i.e., visible light)

        if (fractionalCloudCover < 0.3)
            par0noon = irrAtNoonSurfaceEarthNoCloud*(1-SURFACE_ALBEDO)*(1-iceFrac(iJulianDay))*PAR_FRAC; % W m-2
        else
            par0noon = irrAtNoonSurfaceEarthCloud*(1-SURFACE_ALBEDO)*(1-iceFrac(iJulianDay))*PAR_FRAC; % W m-2
        end
        
        % .................................................................

        % Diurnal surface PAR (Platt et al. 1990): use of a sinusoidal 
        % curve to describe diurnal surface PAR. This curve is symmetric 
        % about noon irradiance (irradiance is maximum at noon ==
        % definition of noon)

        % Sunlight duration (photoperiod)
        sunriseHourAngleRad = acos(-(tan(latitudeRad)*tan(sunDeclinAngleRad)));
        sunriseHourAngleDeg = sunriseHourAngleRad*RAD_TO_DEG; 
        nDaylightHours = (2/15)*sunriseHourAngleDeg; % hours
        if (imag(nDaylightHours) ~= 0)
            nDaylightHours = 24;
        end
        
        sunriseTime = 12 - nDaylightHours*0.5; % h, noon time minus half the hours of daylight, hour from midnight (midnight = 00:00)
        sunsetTime = 12 + nDaylightHours*0.5; % h, noon time plus half the hours of daylight

%         figure
%         hourInDay = [0:1:24];
%         PARsurfaceOcean = PAR0noon * sin(pi*(hourInDay-sunriseTime)/nDaylightHours); % W m-2, diurnal surface PAR    
%         plot(hourInDay,PARsurfaceOcean)
        
        % Daily surface PAR
        par0t = @(t,par0noon,sunriseTime,nDaylightHours) par0noon*sin(pi*(t-sunriseTime)/nDaylightHours); % W m-2
        
        % .................................................................
        
        % Evaluate the integral from sunriseTime to sunsetTime at the 
        % given PAR0noon, sunriseTime and nDaylightHours. This way, we get 
        % the total amount of radiation per m2 in a day. We then average 
        % this quantity per hour and then per second. The integral in time 
        % of W m-2 (=J s-1 m-2) is ([J s-1 m-2] x [s]) J m-2
        totPar0photoperiod = 3600*integral(@(t)par0t(t,par0noon,sunriseTime,nDaylightHours), sunriseTime, sunsetTime); % J m-2
        
        % totPAR0photoperiod is the total amount of sunlight received per unit area in
        % time sunsetTime-sunriseTime (daylight period)

        % PAR0photoperiod is the average amount of sunlight received during the
        % daylight period.
        par0photoperiod(iJulianDay) = totPar0photoperiod/(nDaylightHours*3600); % J s-1 m-2 = W m-2
        lengthPhotoperiod(iJulianDay) = nDaylightHours;
        
        
    else
        
        par0photoperiod(iJulianDay) = 0; % W m-2
        lengthPhotoperiod(iJulianDay) = 0; % h
        
    end
    
end % end loop through iJulianDay

end % function calculatePAR0
