# MATLAB Tools to Work with Oceanographic Datasets

This repository contains a collection of MATLAB scripts that I have developed for reading, formatting and visualising freely-available oceanographic datasets stored in `.nc` (NetCDF) files. It includes tools to organise data into a four-dimensional structure: `latitudes x longitudes x depth levels x 12 months`, essentially creating a climatology. This structure facilitates further data processing, such as data-model validation or model input preparation. Developed across various projects, I actively maintain and update these scripts to meet evolving research needs.

![README_cover](https://github.com/user-attachments/assets/44fb8622-59b6-4843-b2db-d3a930cdb445)

## Requirements

To use the content of this repository, ensure you have the following.
- [MATLAB](https://mathworks.com/products/matlab.html) version R2021a or later installed.
- Third-party functions downloaded from [MATLAB's File Exchange](https://mathworks.com/matlabcentral/fileexchange/): `worstcase`, `m_map`, `brewermap` and `subaxis`. Once downloaded, place these in the `./resources/external/` directory.
- [CO2SYS algorithm](https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/CO2SYS/co2rprt.html) for carbonate system calculations. Once downloaded, place it in the `./resources/external/` directory.
- `copernicusmarine` toolbox, essential for downloading oceanographic data from the Copernicus Marine Environmental Monitoring Service (CMEMS). Installation instructions are in this related [repository](https://github.com/annarufas/jupyter-matlab-ocean-satellite-data-toolbox).
- Jupyter Notebook, included in the latest [Anaconda Python distribution](https://www.anaconda.com/).

## Repository Structure

 - `./code/`: contains the MATLAB scripts (and additional Jupyter Notebooks and shell scripts) for downloading, reading, formating and visualising data (*provided, see "Scripts Overview"*).
 - `./data/`
    - `./raw/`: raw data downloaded from source URLs (*not provided, see below for details*).
    - `./processed/`: processed data generated by MATLAB scripts (*not provided, see below for details*).
- `./resources/`
    - `./external/`: third-party resources for plotting and functions (*see "Requirements" section*).
    - `./internal/`: custom MATLAB functions generated specifically for plotting (*provided*).
- `./figures/`: figures generated from processed data (*provided*).

Due to large file sizes and variety of licenses that limit re-distribution of data in various ways, raw data are not hosted in the `./data/raw/` folder. Instead, the links for manually obtaining these data (as `.nc` files) are provided in the "Data Sources" section below as well as within the MATLAB scripts in the `./code/` folder. If manual access is difficult, we iclude scripts to download data programmatically. Note that data URLs may change over time, potentially interrupting access. The processed data created by the MATLAB scripts (`.mat`) are placed in the `./data/processed/` folder but cannot be provided due to their large size.

## Data Sources

The oceanographic datasets used in this repository are sourced from the following open-access resources:

| Open-access resource                                | Sector                             |
|-----------------------------------------------------|------------------------------------|
| Copernicus Marine Environmental Monitoring Service ([CMEMS](https://data.marine.copernicus.eu/products)) Data Store | Ocean biogeochemistry, ocean physics and ocean colour data |
| Earth System Grid Federation ([ESGF](https://aims2.llnl.gov/search/cmip6/?mip_era=CMIP6)) | Climate data from Coupled Model Intercomparison Projects phase 6 (CMIP6) |
| ESA Biological Pump and Carbon Exchange Processes ([BICEP](https://bicep-project.org/About.html)) project | Ocean colour and ocean  biogeochemistry data |
| ESA Ocean Colour Climate Change Initiative ([OC-CCI](https://www.oceancolour.org)) | Ocean colour data |
| General Bathymetric Chart of the Ocean ([GEBCO](https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global)) | Bathymetry data |
| Global Ocean Data Analysis Project ([GLODAP](https://glodap.info)) | Ocean biogeochemistry, physics and carbonate system variables |
| Institut Français de Recherche pour l'Exploitation de la Mer (IFREMER) [Mixed Layer Depth Climatology](https://cerweb.ifremer.fr/deboyer/mld/home.php) website | Mixed layer depth data |
| NASA [Ocean Color](https://oceancolor.gsfc.nasa.gov) website | Ocean colour data |
| NOAA National Centers for Environmental Information (NCEI) [AVHRR Pathfinder SST](https://www.ncei.noaa.gov/products/avhrr-pathfinder-sst) | Sea surface temperature data |
| NOAA National Centers for Environmental Information (NCEI) [World Ocean Atlas](https://www.ncei.noaa.gov/products/world-ocean-atlas) | Ocean biogeochemistry and ocean physics data |
| Oregon State University's [Ocean Productivity Site](http://orca.science.oregonstate.edu/npp_products.php) | Net primary production data |

## Datasets Used

The specific oceanographic datasets accessed and processed by this repository include:

- Aeolian dust deposition
    - CMIP6 [NCAR-CESM2 historical simulation](https://aims2.llnl.gov/search/cmip6/)
- Bathymetry 
    - [GEBCO](https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global)
- Carbonate system variables (carbonate ion concentration, omega calcite and omega aragonite)
    - [GLODAPv2.2016b](https://www.nodc.noaa.gov/archive/arc0107/0162565/2.2/data/0-data/mapped/)/[CO2SYS](https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/oceans/CO2SYS/co2rprt.html)
- Chlorophyll a (chla) concentration
    - [NASA Aqua-MODIS sensor](https://oceancolor.gsfc.nasa.gov/about/missions/aqua/)
    - [OC-CCI](https://www.oceancolour.org/thredds/ncss/cci/v6.0-release/geographic/monthly/chlor_a/)
    - [CMEMS global BGC reanalysis](https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description)
- Diffuse attenuation coefficient (k<sub>d</sub>)
     - [NASA Aqua-MODIS sensor](https://oceancolor.gsfc.nasa.gov/about/missions/aqua/)
- Euphotic layer depth (z<sub>eu</sub>)
    - Calculated from k<sub>d</sub> from [CMEMS global BGC reanalysis](https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description) and [NASA Aqua-MODIS sensor](https://oceancolor.gsfc.nasa.gov/about/missions/aqua/)
- Mixed layer depth (MLD)
    - [CMEMS global PHYS reanalysis](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description)
    - [IFREMER](https://cerweb.ifremer.fr/deboyer/data/mld_DReqDTm02_c1m_reg2.0.nc)
- Net primary production (NPP)
    - [Oregon State University's Ocean Productivity Site](http://orca.science.oregonstate.edu/npp_products.php) (VGPM, CbPM, CAFE)
    - [BICEP](https://catalogue.ceda.ac.uk/uuid/69b2c9c6c4714517ba10dab3515e4ee6/)
    - [CMEMS global BGC reanalysis](https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description)
- Nutrients (nitrate, silicate, phosphate and dissolved oxygen concentration)
    - [World Ocean Atlas 2023](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/)
- Photosynthetic active radiation at the surface ocean (PAR<sub>0</sub>)
    - [NASA Aqua-MODIS sensor](https://oceancolor.gsfc.nasa.gov/about/missions/aqua/)
    - [NASA SeaWiFS sensor](https://oceancolor.gsfc.nasa.gov/about/missions/seawifs/)
    - Calculated using astronomic/trigonometric equations and data inputs of sea ice fraction from [CMEMS global PHYS reanalysis](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) and cloud cover fraction from [Pincus et al. (2008)](https://doi.org/10.1029/2007JD009334)
- Sea surface temperature (SST)
    - NOAA National Centers for Environmental Information (NCEI) [AVHRR Pathfinder SST](https://www.ncei.noaa.gov/products/avhrr-pathfinder-sst)
- Seawater temperature
    - [World Ocean Atlas 2023](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/)
    - [CMEMS global PHYS reanalysis](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description)
- Mesozooplankton concentration
    - CMIP6 [IPSL-PISCES historical simulation](https://aims2.llnl.gov/search/cmip6/)
    - CMIP6 [GFDL-COBALT historical simulation](https://aims2.llnl.gov/search/cmip6/)
    - CMIP6 [UKESM-MEDUSA historical simulation](https://aims2.llnl.gov/search/cmip6/)

## Scripts Overview

The following scripts are available in the `./code/` folder. 

| Num| Script name                                    | Script action                                           |
|----|------------------------------------------------|----------------------------------------------------------
| 1  | downloadBGCandPHYSfromCMEMS.ipynb              | Downloads data from CMEMS (must be run before script 5)  |                     
| 2  | downloadChlaFromOCCCI.m                        | Downloads data from OC-CCI (must be run before script 7) |  
| 3  | ncreadAerosolDustDepositionFromCMIP6.m         | Creates `dustflux_cmip6_ncarcesm2.mat` (192 x 288 x 12) |
| 4  | ncreadBathymetryFromGEBCO.m                    | Creates `bathymetry_gebco.mat` (1080 x 2160)            |
| 5  | ncreadBGCandPHYSfromCMEMS.m                    | Creates `chla_cmems_bgc.mat`, `kd_cmems_bgc.mat`, `mld_cmems_phys.mat`, `icefrac_cmems_phys.mat` (1080 x 2160 x 12) and `temp_cmems_phys.mat` (1080 x 2160 x 50 x 12) |
| 6  | ncreadBGCandPHYSfromWOA.m                      | Creates `nit_monthly_woa23.mat`, `phos_monthly_woa23.mat`, `sil_monthly_woa23.mat`, `temp_monthly_woa23.mat`, `oxy_monthly_woa23.mat` (180 x 360 x 102 x 12) and `temp_annual_woa23.mat` (180 x 360 x 102) |
| 7  | ncreadChlaFromNASAandOCCCI.m                   | Creates `chla_aquamodis.mat` (4320 x 8640 x 12), `chla_seawifs.mat` (2160 x 4320 x 12) and `chla_occci.mat` (4320 x 8640 x 12) | 
| 8  | ncreadCloudCoverFromPincus.m                   | Creates `cloudcover_pincus.mat` (72 x 144 x 12)         |
| 9  | ncreadKdFromNASA.m                             | Creates `kd_aquamodis.mat` (4320 x 8640 x 12)           |
| 10 | ncreadNPPfromBICEP.m                           | Creates `npp_bicep.mat` (2160 x 4320 x 12)              |
| 11 | ncreadNPPfromOceanProductivitySite.m           | Creates `npp_cafe_seawifs.mat`, `npp_cafe_modis.mat`, `npp_cbpm_modis.mat` and `npp_vgpm_modis.mat` (1080 x 2160 x 12)  |
| 12 | ncreadMLDfromIFREMER.m                         | Creates `mld_ifremer.mat` (90 x 180 x 12)               |
| 13 | ncreadPAR0fromNASA.m                           | Creates `par0_aquamodis.mat` (4320 x 8640 x 12) and `par0_seawifs.mat` (2160 x 4320 x 12)  |
| 14 | ncreadSSTfromPathfinder.m                      | Creates `sst_pathfinder_v5.mat` (4096 x 8192 x 12)      |
| 15 | ncreadZooplanktonFromCMIP6.m                   | Creates `mesozoo_cmip6_pisces.mat` (64 x 128 x 75), `mesozoo_cmip6_cobalt.mat` (180 x 360 x 35) and `mesozoo_cmip6_medusa.mat` (64 x 128 x 75) |
| 16 | createGridFromBathymetricData.m                | Creates `grid_GEBCO_2160_1080.mat` (1080 x 2160 x 500) and `grid_GEBCO_360_180.mat` (180 x 360 x 500) (run after script 4) |
| 17 | createGriddedCarbonateSystemClimatology.m      | Creates `co3ion_co2sys.mat`, `omegacalcite_co2sys.mat` and `omegaaragonite_co2sys.mat` (180 x 360 x 33 x 12) |
| 18 | createGriddedPAR0climatology.m                 | Creates `par0_monthly_calculated.mat` (180 x 360 x 12) and `par0_daily_calculated.mat` (180 x 360 x 365) (run after scripts 6, 9 and 16) |
| 19 | createGriddedZeuClimatology.m                  | Creates `zeu_calculated_kdcmems_mldcmems_pointonepercentpar0.mat`, `zeu_calculated_kdaquamodis_mldcmems_pointonepercentpar0.mat` (1080 x 2160 x 12) (run after scripts 5 and 9) |
| 20 | createGriddedNPPclimatologyFromCarrAlgorithm.m | Creates `npp_carr2002_seawifs_pathfinder.mat` (180 x 360 x 12) |
| 21 | regridZooplanktonConcentrationFromCMIP6.m      | Called by script 15                                 |
| 22 | calculatePAR0fromTrigonometricEquations.m      | Called by script 18                                 |
| 23 | Carr2002algorithm.m                            | Called by script 20                                 |
| 24 | processNASAsensorData.m                        | Called by scripts 8 and 13                          |
| 25 | prepareDataForPlotting.m                       | Creates figures to show monthly climatological data (figures with `_monthly_` infix)|
| 26 | plotCrossSourceComparisonMaps.m                | Creates figures to show comparisons of the same variable across datasets (figures with `_comparison_` infix) |
| 27 | compareInterpolationMethods.m                  | Applies two different interpolation methods to fill data gaps (especially relevant in polar latitudes) and visualises the output |
| 28 | submit_zoo_regridding.sh                       | Submits script 21 to the SLURM job scheduler |   

