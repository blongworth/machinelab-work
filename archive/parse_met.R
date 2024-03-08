#!/usr/bin/R
# concatenate LECS SD data files

library(tidyverse)
library(mlabtools)
library(scattermore)

# list files
file_dir <- "C:/Users/brett/Desktop/LECS_2023/LECS_Surface/SDHC"
files <- list.files(file_dir, pattern = "^2023", full.names = TRUE)

met <- parse_met(files)

write_csv(status, "surface_met.csv")

ggplot(met, aes(ts, par)) +
  geom_scattermore()

met |> 
  filter(ts > "2023-05-08") |> 
ggplot(aes(ts, par)) +
  geom_scattermore()

ggplot(met, aes(ts, wind_speed)) +
  geom_scattermore()

met |> 
  filter(ts > "2023-05-08") |> 
ggplot(aes(ts, wind_speed)) +
  geom_scattermore()

ggplot(met, aes(ts, wind_dir)) +
  geom_scattermore()
