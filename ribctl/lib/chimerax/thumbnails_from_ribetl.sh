#!/bin/bash

# Check if data directory path is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 /path/to/RIBETL_DATA"
    exit 1
fi

RIBETL_DATA="$1"

# Maximum quality GIF conversion with compatible settings
convert_to_gif() {
    local input=$1
    local output=$2
    
    echo "Converting $input to $output with maximum quality settings"
    
    ffmpeg -y -i "$input" -filter_complex "\
        fps=30,\
        scale=720:-1:flags=bicubic,\
        crop=iw*0.5:ih:iw*0.25:0,\
        split[s0][s1];\
        [s0]palettegen=max_colors=256:reserve_transparent=0:stats_mode=full[p];\
        [s1][p]paletteuse=dither=none\
    " -loop 0 "$output"
}

# Extract a single frame with specific quality and resolution settings
extract_frame() {
    local input=$1
    local output=$2
    local timestamp=${3:-0}  # Default to first frame if no timestamp provided
    
    echo "Extracting frame at timestamp $timestamp from $input with matched quality settings"
    
    ffmpeg -y -i "$input" \
        -ss "$timestamp" \
        -vframes 1 \
        -vf "scale=720:-1:flags=bicubic,crop=iw*0.5:ih:iw*0.25:0" \
        -compression_level 0 \
        "$output"
}

# Function to process a single MP4 file
process_video() {
    local input="$1"
    local gif_output="${input%.*}.gif"
    local png_output="${input%.*}.png"
    
    # Only process if outputs don't exist
    if [ ! -f "$gif_output" ]; then
        convert_to_gif "$input" "$gif_output"
    fi
    
    if [ ! -f "$png_output" ]; then
        extract_frame "$input" "$png_output"
    fi
}

# Export functions for GNU parallel
export -f convert_to_gif
export -f extract_frame
export -f process_video

# Find all MP4 files and process them in parallel
find "$RIBETL_DATA" -type f -name "*.mp4" | \
    parallel --bar --jobs 75% process_video

echo "All conversions completed!"