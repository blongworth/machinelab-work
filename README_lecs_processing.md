---
title: LECS data processing
author: Brett Longworth
date: 2024-03-18
---

# Raw files

Raw lecs data is recorded on the surface teensy SD card
in files that contain 4h of data.
These files contain several line types,
defined by the start of the line.
Lines starting with `D:` contain ADV data lines, 
with velocity and other ADV data, and rinko temp and O2.
`S:` lines contain ADV status data, including a timestamp
from the lander teensy RTC and one from the ADV RTC.
`M:` lines contain surface met data with a timestamp.

# Initial processing

Raw lecs data is read and processed by files
before joining files into a larger dataset.
Lines are given a row index then divided into their types
before parsing into tabular data.
Rough calibrations are applied to various data, 
including temp, O2, pH, pressure, 
ADV parameters like voltage and inertial measurements, etc.
ADV data is timestamped using the timestamp of the closest status line,
using ADV count to account for missing packets.
There are cases where it is impossible to give the correct timestamp.
When possible, these are assigned NA, but in the case of skipped packets,
the byte integer count makes it impossible to know how many counts are skipped.
The data are filtered to remove rows with impossible values
in initial processing.

# Initial filtering

Before calibration, the complete dataset is filtered
to remove unrealistic values in things like pH, temp,
and depth.

# Calibration

Data have an initial calibration applied using empirical
and sensor constants along with data from lab calibrations.
Fine calibration and drift correction uses dataset
from a co-deployed SeapHOx, deployed roughly monthly
for a week at a time. 
A linear model is generated for data from the co-deployment periods.
pH regresses lecs raw pH and temp onto temperature corrected
SeapHOx pH.
To correct sensor drift, the fit coefficients for a given time are
fitted linearly between calibrations. 

# Aggregations

For easier long-term trend plotting, hourly and minute-aggregated data
are stored.

# Fluxes

Fluxes of heat, H+, and O2 are calculated by determining the 
integrated cross power spectral density between a parameter 
and the vertical velocity for the time aggregation of interest
(hourly). Spectral data are filtered to remove longer period
fluxes due to wave and tidal movement,
and shorter periods that are below the nyquist frequency
or just noise.

# files

* Raw data are stored in `data/SD Card Data`
* Processed data are stored in `data/processed`
* Initally processed data are in Parquet format in files ending in `all`
* Filtered, calibrated data are in Parquet files ending in `filtered_calibrated`
* Aggregated data are in csv and Parquet files ending in `minute` or `hourly`

# Code 

Processing is done in R using functions from `https://github.com/blongworth/mlabtools`
and scripts from `https://github.com/blongworth/machinelab-work`.
Code by Manish Devana for data processing in python is in
`https://github.com/WHOIGit/LECS_MachineLAB`. 
Spectral methods are from work by Long et al. `https://doi.org/10.1029/2020JC016637`

# Procedure

1. Run `lecs_process_all.R`. This will generate met, status, 
and ADV parquet datasets for each directory of LECS raw files.
This will also filter impossible data and timestamp ADV data

2. Parquet files are written such that Arrow views them as a single dataset,
partitioned by download. 
Combining them explicitly with `combine.qmd` is no longer necessary.

3. Run `filter.qmd` to filter bad data by value. 
Inspection of new data is needed to see what kind of weird data 
is in there.

4. Run `calibrate.qmd` to calibrate pH using seapHOx data.
Add script sections as needed to incorporate new cal data.

5. Run `resample.qmd` to create hourly and minute-averaged data
and store as parquet and csv.

Additional scripts are needed for imputation and handling
of missing data and calculation of hourly fluxes.