#!/usr/bin/R
# concatenate LECS SD data files

library(tidyverse)
library(lubridate)
library(scattermore)


# read files to one df
read_status <- function(file) {
   lines <- readLines(file)
   lines[grep('S:', lines)]
}

parse_status <- function(files) {
  status_lines <- map(files, read_status, .progress = TRUE) |> 
    reduce(c) |> 
    str_remove("^.*S:")
  
  status = read_csv(I(status_lines), col_names = c("hour", "min", "sec", "day", "month", "year", 
                                                   "vmin", "vsec", "vday", "vhour", "vyear", "vmo",
                                                   "bat", "ss", "head", "pitch", "roll", "temp", 
                                                   "empty", "CR", "BV", "PWR")) |> 
    filter(year == 2023) |>
    mutate(ts = make_datetime(year = year, month = month, day = day,
                              hour = hour, min = min, sec = sec),
           bat = bat * .1)
}

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
