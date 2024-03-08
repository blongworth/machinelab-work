## Get and process data

library(mlabtools)

# Make list of all files
file_dir <- "data/SD Card Data/LECS_surface_sd/lecs_surface_2023-11-15"
files <- list.files(file_dir, pattern = "^202[3|4]", full.names = TRUE)

# Process files into a list containing data frames for ADV, status, and Met

plan(multicore)

tic()
lecs_data_p <- lecs_parse_files_p(files)
toc()
tic()
write_csv(lecs_data[["met"]], "lecs_met_2023-11-15.csv")
write_csv(lecs_data[["status"]], "lecs_status_2023-11-15.csv")
write_csv(lecs_data[["adv_data"]], "lecs_adv_2023-11-15.csv")
toc()

# Make list of all files

file_dir <- "data/SD Card Data/LECS_surface_sd/lecs_surface_2023-09-05/"
files <- list.files(file_dir, pattern = "^202[3|4]", full.names = TRUE)

# Process files into a list containing data frames for ADV, status, and Met

plan(multicore)

tic()
lecs_data_p <- lecs_parse_files_p(files)
toc()
tic()
write_csv(lecs_data[["met"]], "lecs_met_2023-09-05.csv")
write_csv(lecs_data[["status"]], "lecs_status_2023-09-05.csv")
write_csv(lecs_data[["adv_data"]], "lecs_adv_2023-09-05.csv")
toc()