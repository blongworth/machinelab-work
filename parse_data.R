#!/usr/bin/R
# concatenate LECS SD data files

library(tidyverse)
library(lubridate)
library(scattermore)


# read files to one df
read_data <- function(file) {
   lines <- readLines(file)
   lines[grep('^D', lines)]
}

parse_data <- function(files) {
  data_lines <- map(files, read_data, .progress = TRUE) |> 
    reduce(c) |> 
    str_remove("^D:")
  
  dl = read_csv(I(data_lines),
                col_names = c("count", "pressure", 
                              "vx", "vy", "vz",
                              "a1", "a2", "a3", 
                              "corr1", "corr2", "corr3",
                              "ana_in", "ana_in2", "pH",
                              "temp", "DO"))
}

# list files
file_dir <- "C:/Users/brett/Desktop/LECS_2023/LECS_Lander/SDHC"
files <- list.files(file_dir, pattern = "^2023", full.names = TRUE)

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
