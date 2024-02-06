#!/usr/local/bin/Rscript
# read and parse all LECS SD card data to csv

library(mlabtools)
library(tictoc)

out_dir <- "data/processed/surface/"

tic()
lecs_process_data("2023-09-05", "/Users/brett/Projects/machinelab-work/data/SD Card Data/LECS_surface_sd/lecs_surface_2023-09-05/clean", out_dir, dedupe = TRUE)
lecs_process_data("2023-11-15", out_dir = out_dir, dedupe = TRUE)
lecs_process_data("2024-01-23", out_dir = out_dir, dedupe = TRUE)
toc()
