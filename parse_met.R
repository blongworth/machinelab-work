#!/usr/bin/R
# concatenate LECS SD data files

library(tidyverse)
library(lubridate)
library(scattermore)


# read files to one df
read_met <- function(file) {
   lines <- readLines(file)
   lines[grep('M:', lines)]
}

parse_met <- function(files) {
  met_lines <- map(files, read_met, .progress = TRUE) |> 
    reduce(c) |> 
    str_remove("^.*M:")
  
  met = read_csv(I(met_lines), col_names = c("hour", "min", "sec", "day", "month", "year", 
                                                   "par", "wind_speed", "wind_dir")) |> 
    filter(year == 2023) |>
    mutate(ts = make_datetime(year = year, month = month, day = day,
                              hour = hour, min = min, sec = sec))
}

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
