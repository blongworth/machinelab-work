#!/bin/bash

# Default range
start=1
end=15

# Check for optional arguments
if [[ $# -eq 2 ]]; then
  start=$1
  end=$2
elif [[ $# -ne 0 ]]; then
  echo "Usage: $0 [start] [end]"
  echo "       start and end are optional and must both be provided to define the range."
  exit 1
fi

# Set the base name for the files
data_prefix="lecs_ml_data_"
timestamp_prefix="lecs_ml_timestamps_"
file_extension=".dat"

# Loop through numbers 1 to 15

for i in $(seq $start $end); do
  # Construct the filenames
  data_file="${data_prefix}${i}${file_extension}"
  timestamp_file="${timestamp_prefix}${i}${file_extension}"

  # Check if both files exist
  if [[ -f "$data_file" && -f "$timestamp_file" ]]; then
    echo "Processing files: $data_file and $timestamp_file"

    # Run MATLAB with the script and pass the file names as arguments
    matlab -nodisplay -nosplash -r "eddyflux_batch('$data_file', '$timestamp_file'); exit;"
  else
    echo "Files $data_file or $timestamp_file not found. Skipping."
  fi
done