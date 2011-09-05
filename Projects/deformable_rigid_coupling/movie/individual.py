from video import *
import sys

name=["poisson-300","poisson-400","poisson-450","poisson-490","poisson-499","poisson-500",
      "stiff_0.5","stiff_3","stiff_20","stair-good","armadillo2","armadillo5",
      "armadillo-stairs-stiffness-2","armadillo-stairs-stiffness-4","armadillo-stairs-stiffness-8","tori/hires"]
    
for i in name:
    video_dir="cloth_video_tmp"
    if os.path.isdir(video_dir):
        shutil.rmtree(video_dir)

    video=VIDEO(video_dir)
    title_length=video.fps*5*1.2
    caption_length=int(video.fps*3.5*1.2)
    caption_long_length=int(video.fps*5*1.2)
    caption_short_length=int(video.fps*2*1.2)
    caption_tiny_length=int(video.fps*1*1.2)
    video.add_directory(i)
    video.make_movie(i+'-movie')

