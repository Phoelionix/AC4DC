#!/bin/zsh
mkdir -p video
rm video/*
python3 scripts/plot_slice_video.py Carbon_HR1 600 50
ffmpeg -r 30 -f image2 -s 1920x1080 -i video/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p out.mp4