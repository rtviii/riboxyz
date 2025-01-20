#!/bin/bash

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

# Usage examples:
# Extract first frame:
# ./script.sh input.mp4 output.png
#
# Extract frame at specific timestamp:
# ./script.sh input.mp4 output.png 1.5  # Extract frame at 1.5 seconds

# Quality settings matched from GIF converter:
# - 720p width with maintained aspect ratio
# - Bicubic scaling for high quality
# - Center cropped to middle 50% of frame
# - Maximum quality PNG output (compression_level 0)

# Execute the extraction if arguments are provided
if [ $# -ge 2 ]; then
    extract_frame "$1" "$2" "${3:-0}"
else
    echo "Usage: $0 input_file output_file [timestamp]"
    echo "Example: $0 video.mp4 frame.png 1.5"
    exit 1
fi