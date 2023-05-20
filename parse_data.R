#!/usr/bin/R
# concatenate LECS SD data files

library(tidyverse)
library(scattermore)
library(mlabtools)

# list files
file_dir <- "C:/Users/brett/Desktop/LECS_2023/LECS_Lander/SDHC"
files <- list.files(file_dir, pattern = "^2023", full.names = TRUE)

# read files to one df
df <- parse_data(files[541])

write_csv(df, "lander_data.csv")

ggplot(df, aes(1:nrow(df), pH)) +
  geom_scattermore()

ggplot(df, aes(1:nrow(df), temp)) +
  geom_scattermore()

ggplot(df, aes(1:nrow(df), DO)) +
  geom_scattermore()

ggplot(df, aes(1:nrow(df), corr1)) +
  geom_scattermore()

