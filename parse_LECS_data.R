#!/usr/bin/R
# concatenate LECS SD data files

library(tidyverse)
library(mlabtools)

# Lander data
# list files
file_dir <- "data/SD Card Data/LECS_surface_sd/lecs_surface_2024-01-23/"
files <- list.files(file_dir, pattern = "^2023", full.names = TRUE)


status <- parse_status(files)
write_csv(status, "data/LECS_2023/lander_status.csv")

# maybe need to do data in chunks?
adv_data <- parse_data(files)
write_csv(adv_data, "data/LECS_2023/lander_data.csv")


# Surface SD data
# list files
file_dir <- "data/LECS_2023/LECS_Surface/SDHC/"
files <- list.files(file_dir, pattern = "^2023", full.names = TRUE)

# status
status <- parse_status(files)
write_csv(status, "data/LECS_2023/surface_status.csv")

# ADV data
# maybe need to do data in chunks?
adv_data <- parse_data(files[130])
write_csv(adv_data, "data/LECS_2023/surface_data.csv")

# met data
met <- parse_met(files)
write_csv(met, "data/LECS_2023/surface_met.csv")
