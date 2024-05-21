#!/usr/local/bin/Rscript
# read and parse all LECS SD card data to csv

library(mlabtools)
library(tictoc)

out_dir <- "data/processed/surface/"

tic()
lecs_process_data(date = "2023-09-05", file_dir = "/Users/brett/Projects/machinelab-work/data/SD Card Data/LECS_surface_sd/lecs_surface_2023-09-05/clean", out_dir = out_dir, dedupe = FALSE)
lecs_process_data(date = "2023-11-15", file_dir = "/Users/brett/Projects/machinelab-work/data/SD Card Data/LECS_surface_sd/lecs_surface_2023-11-15/clean", out_dir = out_dir, dedupe = FALSE)
lecs_process_data(date = "2024-01-23", file_dir = "/Users/brett/Projects/machinelab-work/data/SD Card Data/LECS_surface_sd/lecs_surface_2024-01-23/clean", out_dir = out_dir, dedupe = FALSE)
lecs_process_data(date = "2024-02-28", file_dir = "/Users/brett/Projects/machinelab-work/data/SD Card Data/LECS_surface_sd/lecs_surface_2024-02-28/clean", out_dir = out_dir, dedupe = FALSE)
lecs_process_data(date = "2024-05-14", file_dir = "/Users/brett/Projects/machinelab-work/data/SD Card Data/LECS_surface_sd/lecs_surface_2024-05-14/clean", out_dir = out_dir, dedupe = FALSE)
#lecs_process_data("2023-09-05", out_dir = out_dir, dedupe = TRUE)
#lecs_process_data("2023-11-15", out_dir = out_dir, dedupe = FALSE)
#lecs_process_data("2024-01-23", out_dir = out_dir, dedupe = FALSE)
toc()

out_dir <- "data/processed/lander/"

#tic()
#lecs_process_data(date = "2023-07-24", "/Users/brett/Projects/machinelab-work/data/SD Card Data/LECS_lander_sd/lecs_lander_2023-07-24/clean", out_dir = out_dir, dedupe = FALSE)
#toc()
