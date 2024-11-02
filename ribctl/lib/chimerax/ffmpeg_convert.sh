echo "Attempting to convert $1 to $2"
ffmpeg -y -i $1 -filter_complex "fps=20,scale=360:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=128[p];[s1][p]paletteuse=dither=bayer" -loop 0 $2