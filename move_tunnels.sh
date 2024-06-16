#!/bin/bash

# Set the source and destination root directories
SRC_ROOT="/home/rtviii/dev/riboxyz/ribctl/assets/exit_tunnel_work"
DEST_ROOT="/home/rtviii/dev/RIBETL_DATA"

# Ensure the source directory exists
if [ ! -d "$SRC_ROOT" ]; then
    echo "Error: Source directory does not exist: $SRC_ROOT"
    exit 1
fi

# Ensure the destination directory exists
if [ ! -d "$DEST_ROOT" ]; then
    echo "Error: Destination directory does not exist: $DEST_ROOT"
    exit 1
fi

# Loop through all directories in the source
for src_dir in "$SRC_ROOT"/*/; do
    # Extract the 4-character ID
    rcsb_id=$(basename "$src_dir")

    # Check if it's a 4-character directory
    if [[ ! $rcsb_id =~ ^[A-Za-z0-9]{4}$ ]]; then
        echo "Skipping non-standard directory: $rcsb_id"
        continue
    fi

    # Create the destination directory if it doesn't exist
    dest_dir="$DEST_ROOT/$rcsb_id"
    if [ ! -d "$dest_dir" ]; then
        mkdir -p "$dest_dir"
        echo "Created destination directory: $dest_dir"
    fi

    # Create the TUNNELS subdirectory
    tunnels_dir="$dest_dir/TUNNELS"
    mkdir -p "$tunnels_dir"

    # Copy the contents
    if cp -R "$src_dir"* "$tunnels_dir"; then
        echo "Copied contents from $src_dir to $tunnels_dir"
    else
        echo "Error copying from $src_dir to $tunnels_dir"
    fi
done

echo "Operation completed."
