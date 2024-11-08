function [npp] = Carr2002algorithm(chla, sst, par0, isZeuCarr)

% CARR2002ALGORITHM Computes net primary production using the Carr (2002) 
% algorithm (https://doi.org/10.1016/S0967-0645(01)00094-7) for a 
% lat x lon x 12 array, accounting for chlorophyll a concentration, 
% sea surface temperature and photosynthetically active radiation.
% 
%   INPUT: 
%       chla      - chlorophyll a concentration (mg m-3), lat x lon x 12 array
%       sst       - sea surface temperature (ºC), lat x lon x 12 array
%       par0      - photosynthetic active radiation surface ocean (W m-2), lat x lon x 12 array
%       isZeuCarr - Boolean indicating choice of euphotic layer depth (zeu)
%                   calculation (1=Carr 2002, 2=Behrenfeld & Falkowski, 1997)
%
%   OUTPUT:
%       npp - net primary production (mg C m-2 d-1), lat x lon x 12 array
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 4 Nov 2024 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Calculate light-related properties

% Initialise output arrays
[nr_chla,nc_chla,nt_chla] = size(chla);
kPar    = NaN(nr_chla,nc_chla,nt_chla); % attenuation coefficient for PAR 
parZeu  = NaN(nr_chla,nc_chla,nt_chla); % PAR at the euphotic layer depth 
zeuCarr = NaN(nr_chla,nc_chla,nt_chla); % euphotic layer depth, Carr (2002) 
zeuBF   = NaN(nr_chla,nc_chla,nt_chla); % euphotic layer depth,  Behrenfeld & Falkowski (1997) 

% Loop over lat, lon, time
for iRow = 1:nr_chla
    for iCol = 1:nc_chla
        for iMonth = 1:12
            if (~isnan(chla(iRow,iCol,iMonth)))

                % Calculate kPar (Eq. 1 in Carr 2002, after Nelson & Smith, 1991)
                kPar(iRow,iCol,iMonth) = 0.04 + (0.0088.*chla(iRow,iCol,iMonth))... 
                    + 0.054.*chla(iRow,iCol,iMonth).^(0.66); % m-1 

                % Calculate zeu (Eq. 4 in Carr 2002)
                zeuCarr(iRow,iCol,iMonth) = -log(0.01)./kPar(iRow,iCol,iMonth); % m

                % Calculate zeu as in Behrenfeld & Falkowski (1997).
                % Henson et al. (2012) use the Behrenfeld & Falkowski (1997) 
                % algorithm (VGPM model, https://sites.science.oregonstate.edu/ocean.productivity/vgpm.code.php)
                % to compute zeu instead of Carr (2002) proposed equation
                if (chla(iRow,iCol,iMonth) < 1)
                    Ctot=38*(chla(iRow,iCol,iMonth).^0.425); % integrated water column chl
                else
                    Ctot=40.2*(chla(iRow,iCol,iMonth).^0.507);
                end
                zeu = 568.2*(Ctot.^-0.746);
                if (zeu > 102)
                    zeu = 200*(Ctot.^-0.293);
                end
                zeuBF(iRow,iCol,iMonth) = zeu; % m

                % Select euphotic layer depth based on input parameter
                if (isZeuCarr)
                    zeu = zeuCarr;
                else
                    zeu = zeuBF;
                end

                % Calculate PAR at zeu (Eq. 2 in Carr 2002, after Riley 1957)
                parZeu(iRow,iCol,iMonth) = par0(iRow,iCol,iMonth)...
                    .*(1-exp(-kPar(iRow,iCol,iMonth).*zeu(iRow,iCol,iMonth)))...
                    ./(kPar(iRow,iCol,iMonth).*zeu(iRow,iCol,iMonth)); % W m-2 

            end
        end
    end
end

%% Calculate NPP (Eq. 3 in Carr 2002)

ALPHA = 2.64; % mg C (mg chla)-1 d-1 (W m-2)-1
Pmax = 24.*exp(0.09.*sst); % mg C (mg chla)-1 d-1, the Eppley curve (Eppley 1972), after Platt et al. (1991) (see here: https://www.sciencedirect.com/science/article/pii/S0964274997800182?via=ihub)
npp = chla.*(((Pmax.*parZeu)./((Pmax./ALPHA)+parZeu))).*zeu; % mg C m-2 d-1

% The above was coded by S. Henson. I originally coded NPP calculation
% differently (see below). Both codes are equivalent, the difference is in 
% the way they present the Epplye curve (max. phytoplankton growth rate as a 
% function of temperature). Eppley (1972) originally formulated his 
% equation for max. growth rate in units of d-1 (umax). Later on, in 
% 1991, Platt et al. revisited Eppley (1972) equation so that it could 
% be expressed in units of mg C (mg chla)-1 d-1 (Pmax). Platt et al. 
% (1991) equation implicitly assumes a C/chla ratio of (=24/0.59) 
% 40.1 mg C (mg chl)-1. In my calculations, I had originally assumed a 
% C/chla of 50 mg C (mg chl)-1 (after Fasham 1990). 

% % As originally coded by me (umax)
% THETA_C_TO_CHL = 40.1; % 50 g C (g chla)-1 = 0.020 g chla (g C)-1 (Fasham et al. 1990)
% muMax = 0.59.*exp(0.0633.*sst(:,:,:)); % d-1, the Eppley curve (Eppley 1972), after Kremer et al. (2017) – 0.59 d-1 is the max. growth rate at 0ºC
% % muMax = 0.59.*1.88.^(sst(:,:,:)./10) % d-1, the Eppley curve (Eppley 1972), after Kremer et al. (2017) - Q10 factor equivalent of the exponential model
% Ik = muMax(:,:,:)./(2.64/THETA_C_TO_CHL); % W m-2
% chlaProd = chla(:,:,:).*((muMax(:,:,:).*parZeu(:,:,:))./(Ik(:,:,:)+parZeu(:,:,:))); % mg chla m-3 d-1
% npp = chlaProd(:,:,:).*THETA_C_TO_CHL.*zeu(:,:,:); % mg chla m-3 d-1 --> mg C m-2 d-1

%% Optional: visual comparison of euphotic depths

% nSubplots = 2;
% mytitlestring = {'Carr (2002) algorithm','B&F algorithm'};
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.29 0.49],'Color','w')
% haxis = zeros(nSubplots,1);
% 
% for iSubplot = 1:nSubplots  
%     haxis(iSubplot) = subplot(2,1,iSubplot);
%     ax(iSubplot).pos = get(haxis(iSubplot),'Position');
%     if (iSubplot == 1)
%         mydata = zeuCarr;
%     elseif (iSubplot == 2)
%         mydata = zeuBF;
%     end
%     pcolor(flipud(rot90(mydata(:,:,4))))
%     caxis([40 110])
%     shading interp; colormap(jet)
%     set(gca,'xticklabels',[],'yticklabels',[])
%     box on
%     title(sprintf('%s',mytitlestring{iSubplot})) 
%     set(haxis(iSubplot),'FontSize',12)
% end
% ax(1).pos(1) = ax(1).pos(1) - 0.12; ax(1).pos(4) = ax(1).pos(4) + 0.040;
% ax(2).pos(1) = ax(2).pos(1) - 0.12; ax(2).pos(4) = ax(2).pos(4) + 0.040; 
% for iSubplot = 1:nSubplots
%     set(haxis(iSubplot), 'Position', ax(iSubplot).pos) 
% end    
% cb = colorbar('Location','eastoutside');
% cb.Position(1) = 0.85;
% cb.Position(2) = 0.25;
% cb.Position(4) = 0.60; 
% cb.Label.String = 'z_{eu} (m)'; 
% cb.FontSize = 12;

%% Optional: compare with Fig. 4a in Henson et al. (2012)

% nppAnnualMean(:,:) = (365/1000).*mean(npp,3,'omitnan'); % g C m-2 yr-1
% figure(); pcolor(flipud(rot90(nppAnnualMean))); 
% caxis([0 500]); 
% cb = colorbar('FontSize', 15, 'FontWeight', 'bold'); 
% cb.Label.String = 'NPP (g C m^{-2} yr^{-1})';
% shading interp; colormap(jet); box on;

end % Carr2002algorithm