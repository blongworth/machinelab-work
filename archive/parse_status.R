#!/usr/bin/R
# concatenate LECS SD data files

library(tidyverse)
library(mlabtools)
library(scattermore)

# list files
file_dir <- "C:/Users/brett/Desktop/LECS_2023/LECS_Lander/SDHC"
files <- list.files(file_dir, pattern = "^2023", full.names = TRUE)

status <- parse_status(files)

write_csv(status, "lander_status.csv")

status_cl <- status |> 
  filter(bat > 10 & bat < 15)

ggplot(status_cl, aes(ts, bat)) +
  geom_scattermore()

# list files
file_dir <- "C:/Users/brett/Desktop/LECS_2023/LECS_Surface/SDHC"
files <- list.files(file_dir, pattern = "^2023", full.names = TRUE)

status <- parse_status(files)

write_csv(status, "surface_status.csv")

status_cl <- status |> 
  filter(month < 6, 
         bat > 10 & bat < 15)

ggplot(status_cl, aes(ts, bat)) +
  geom_scattermore()
