# cams-act6

Jupyter notebook practicals developed for the 6th joint ECMWF CAMS, EUMETSAT and ESA training in atmospheric composition, held on 16-20 September 2024 at NILU in Kjeller, Norway.

Each directory in this repository contains code, Jupyter notebooks and other (non-data) files related to the practical sessions.

The directories are numbered according to the order in which they appear in [the course programme](https://atmosphere.copernicus.eu/6th-ecmwf-cams-esa-eumetsat-training-atmospheric-composition)

## Preliminary table of contents
This is a preliminary list of the material that might be presented at the meeting training sessions.

1. Satellite: the goal is to show as much data as possible, e.g. S3 Fire Radiative Power, plumes from S3 OCLI, emissions from GOME-2, IASI, TROPOMI NO2 and other species, maybe Scar Burn from S2.
2. In-situ: processing and visualization of in-situ data, e.g. from the Pandora instrument. The instrument is a passive sensor used to monitor the content of some chemical species based on the Differential Optical Absorption Spectroscopy (DOAS) technique.
3. Model: download and visualization of dataset from NWP models, e.g.
   * aerosol optical depth using the dataset [CAMS global atmospheric composition forecasts](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-atmospheric-composition-forecasts?tab=overview)
   * ozone forecast using the dataset [CAMS European air quality forecasts](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-europe-air-quality-forecasts?tab=overview)
4. Emissions: download and map visualization of anthropogenic emissions using the dataset [CAMS global emission inventories](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-emission-inventories?tab=overview)
5. CAMS-MOS: evaluation of climate and air quality models. Comparison of in-situ data with model output statistics (CAMS-MOS) using the MOS dataset 
6. Averaging Kernels: comparison of models with satellite data, e.g. for NO2 from [TROPOMI](https://www.tropomi.eu/), using the averaging kernels statistical algorithm.
7. Aeroval: model evaluation with in-situ data using the AeroVal APIs [CAMS European air quality forecasts optimised at observation sites](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-europe-air-quality-forecasts-optimised-at-observation-sites?tab=overview) and surface measurement e.g. from [AeroVal APIs](https://aeroval.met.no/)
8. Air quality index: calculation of the air quality index forecast from $NO_2$, $O_3$, $SO_2$, $PM10$ and $PM2.5$ concentrations using the dataset [CAMS global atmospheric composition forecasts](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-atmospheric-composition-forecasts?tab=overview)

## Access to the data
By September 16th 2024 The Climate Data Store and the Atmospheric Data Store will be decommissioned and all the requests for data access will be forwarded to the new web services, CDS-Beta and ADS-Beta respectively. In order to access these new web services you need to register at the [ECMWF](https://www.ecmwf.int/). The next step is to register to [CDS-Beta](https://cds-beta.climate.copernicus.eu/) and [ADS-Beta](https://ads-beta.atmosphere.copernicus.eu/). In order to avoid writing the web service url and your access token (key) you may write them in a hidden file **.cdsapirc** in your home folder as follow:
```
url: https://ads.atmosphere.copernicus.eu/api/v2
key: <replace with your access key>
```
