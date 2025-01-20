# echo "Attempting to convert $1 to $2"
# ffmpeg -y -i $1 -filter_complex "fps=20,scale=360:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=128[p];[s1][p]paletteuse=dither=bayer" -loop 0 $2



# echo "Attempting to convert $1 to $2"
# ffmpeg -y -i $1 -filter_complex "fps=20,scale=360:-1:flags=lanczos,crop=iw*0.5:ih:iw*0.25:0,split[s0][s1];[s0]palettegen=max_colors=128[p];[s1][p]paletteuse=dither=bayer" -loop 0 $2

#!/bin/bash

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

# Usage example
# ./script.sh input.mp4 output.gif

# Quality maximizing improvements:
# - 30 FPS for smooth motion
# - 720p width for high resolution
# - Bicubic scaling for quality
# - Center crop to middle 50% of frame
# - Full frame statistics for optimal palette
# - Maximum 256 color palette
# - No dithering for maximum sharpness

# Execute the conversion if arguments are provided
if [ $# -eq 2 ]; then
    convert_to_gif "$1" "$2"
else
    echo "Usage: $0 input_file output_file"
    exit 1
fi