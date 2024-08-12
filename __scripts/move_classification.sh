#!/bin/bash

# Set the source and destination root directories
SRC_DIR="/home/rtviii/dev/riboxyz/ribctl/logs/hmm_classification_reports"
DEST_ROOT="/home/rtviii/dev/RIBETL_DATA"

# Ensure the source directory exists
if [ ! -d "$SRC_DIR" ]; then
    echo "Error: Source directory does not exist: $SRC_DIR"
    exit 1
fi

# Ensure the destination root directory exists
if [ ! -d "$DEST_ROOT" ]; then
    echo "Error: Destination root directory does not exist: $DEST_ROOT"
    exit 1
fi

# Loop through all .json files in the source directory
for src_file in "$SRC_DIR"/*.json; do
    # Check if file exists (in case no .json files are found)
    [ -e "$src_file" ] || continue

    # Extract the RCSB_ID from the filename
    filename=$(basename "$src_file")
    rcsb_id="${filename%.json}"

    # Check if it's a 4-character ID
    if [[ ! $rcsb_id =~ ^[A-Za-z0-9]{4}$ ]]; then
        echo "Skipping file with non-standard ID: $filename"
        continue
    fi

    # Create the new filename
    new_filename="classification_report_${rcsb_id}.json"

    # Set the destination directory
    dest_dir="$DEST_ROOT/$rcsb_id"

    # Check if the destination directory exists
    if [ ! -d "$dest_dir" ]; then
        echo "Destination directory does not exist: $dest_dir"
        continue
    fi

    # Copy and rename the file
    if cp "$src_file" "$dest_dir/$new_filename"; then
        echo "Copied and renamed $filename to $dest_dir/$new_filename"
    else
        echo "Error copying $filename to $dest_dir/$new_filename"
    fi
done

echo "Operation completed."
