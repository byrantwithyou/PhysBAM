#!/usr/bin/python

from video import *
import sys

header_color='#ffaa00'
subtle_color='#808080'

slide("title.png",36,[(36,"Asynchronous Evolution for"),("Implicit Time Integration"),(36,""),(36,"SIGGRAPH 2009 #0189")])

slide("sphere1_explicit.png",24,[(36,"Sphere",header_color),"","Semi-Implicit","(explicit elastic forces)"])
slide("sphere1_implicit.png",24,[(36,"Sphere",header_color),"","Fully-Implicit","","Very damped"])
slide("sphere1_async.png",24,[(36,"Sphere",header_color),"","Asynchronous","","Less damped"])
slide("sphere1_e0.png",36,["e=0.0","",(24,"Little internal support")])
slide("ept25.png",36,["e=0.25"])
slide("ept5.png",36,["e=0.5"])
slide("ept75.png",36,["e=0.75"])
slide("sphere1_e1.png",36,["e=1.0","",(24,"Allows little deformation")])

slide("sphere2_explicit.png",24,[(36,"Two Spheres",header_color),"","Semi-Implicit","(explicit elastic forces)"])
slide("sphere2_implicit.png",24,[(36,"Two Spheres",header_color),"","Fully-Implicit","","poor self-collision handling"])
slide("sphere2_async_rigidpt5.png",24,[(36,"Two Spheres",header_color),"","Asynchronous","","e=0.5"])

slide("armadillo_explicit.png",24,[(36,"Armadillo",header_color),"","Semi-Implicit","(explicit elastic forces)"])
slide("armadillo_implicit.png",24,[(36,"Armadillo",header_color),"","Fully-Implicit","","poor self-collision handling"])
slide("armadillo_async.png",24,[(36,"Armadillo",header_color),"","Asynchronous"])
slide("armadillo_ept5.png",36,["e=0.5","",(24,"Allows deformation")])
slide("armadillo_e1.png",36,["e=1.0","",(24,"Allows little deformation")])

slide("sphere1_bad_explicit.png",24,[(36,"Sphere with",header_color),(36,"a Single Bad Element",header_color),"","Semi-Implicit","(explicit elastic forces)","","Very slow"])
slide("sphere1_bad_implicit.png",24,[(36,"Sphere with",header_color),(36,"a Single Bad Element",header_color),"","Fully-Implicit","","Stable but introduces artifacts"])
slide("sphere1_bad_async_rigidpt5.png",24,[(36,"Sphere with",header_color),(36,"a Single Bad Element",header_color),"","Asynchronous","","e=0.5","","Fast","(8.5x speed-up over semi-implicit)"])

video_dir="async_video_tmp-100"
if os.path.isdir(video_dir):
    shutil.rmtree(video_dir)

video=VIDEO(video_dir)

testing_titles=False
if testing_titles:
    title_length=video.fps*5/10
    caption_length=int(video.fps*3.5/10)
    caption_long_length=int(video.fps*5/10)
    caption_longer_length=int(video.fps*5/10)
    caption_short_length=int(video.fps*2/10)
    caption_tiny_length=int(video.fps*1/10)
    caption_tinier_length=int(video.fps*1/10)
else:
    title_length=video.fps*5*1.2
    caption_length=int(video.fps*3.5*1.2)
    caption_long_length=int(video.fps*5*1.2)
    caption_longer_length=int(video.fps*6*1.2)
    caption_short_length=int(video.fps*2*1.2)
    caption_tiny_length=int(video.fps*1*1.2)
    caption_tinier_length=int(video.fps*1*1.2)
#testing_titles=True
print caption_length

video.add_frame("title.png",title_length)

video.add_frame("sphere1_explicit.png",caption_length)
if not testing_titles: video.add_directory("../render/sphere_1_explicit/Test_24_model_1_explicit_overdamping_1.000000_spheres_1/output")
video.add_frame("sphere1_implicit.png",caption_length)
if not testing_titles: video.add_directory("../render/sphere_1_implicit/Test_24_model_1_implicit_overdamping_1.000000_spheres_1/output",0,120)
video.add_frame("sphere1_async.png",caption_length)
video.add_frame("sphere1_e0.png",caption_length)
if not testing_titles: video.add_directory("../render/sphere_1_asynch/Test_24_model_1_async_overdamping_1.000000_spheres_1_prigid_0.000000/output",0,120)
video.add_frame("ept25.png",caption_short_length)
if not testing_titles: video.add_directory("../render/sphere_1_asynch/Test_24_model_1_async_overdamping_1.000000_spheres_1_prigid_0.250000/output",0,120)
video.add_frame("ept5.png",caption_short_length)
if not testing_titles: video.add_directory("../render/sphere_1_asynch/Test_24_model_1_async_overdamping_1.000000_spheres_1_prigid_0.500000/output",0,120)
video.add_frame("ept75.png",caption_short_length)
if not testing_titles: video.add_directory("../render/sphere_1_asynch/Test_24_model_1_async_overdamping_1.000000_spheres_1_prigid_0.750000/output",0,120)
video.add_frame("sphere1_e1.png",caption_short_length)
if not testing_titles: video.add_directory("../render/sphere_1_asynch/Test_24_model_1_async_overdamping_1.000000_spheres_1_prigid_1.000000/output",0,80)

video.add_frame("sphere2_explicit.png",caption_length)
if not testing_titles: video.add_directory("../render/sphere_2_explicit/Test_24_model_1_explicit_overdamping_1.000000_spheres_2/output",0,150)
video.add_frame("sphere2_implicit.png",caption_length)
if not testing_titles: video.add_directory("../render/sphere_2_implicit/Test_24_model_1_implicit_overdamping_1.000000_spheres_2/output",0,150)
video.add_frame("sphere2_async_rigidpt5.png",caption_length)
if not testing_titles: video.add_directory("../render/sphere_2_asynch/Test_24_model_1_async_overdamping_1.000000_spheres_2_prigid_0.500000/output",0,110)

video.add_frame("armadillo_explicit.png",caption_length)
if not testing_titles: video.add_directory("../render/armadillo_explicit/Test_24_model_2_explicit_overdamping_1.000000/output")
video.add_frame("armadillo_implicit.png",caption_length)
if not testing_titles: video.add_directory("../render/armadillo_implicit/Test_24_model_2_implicit_overdamping_1.000000/output")
video.add_frame("armadillo_async.png",caption_length)
video.add_frame("armadillo_e1.png",caption_short_length)
if not testing_titles: video.add_directory("../render/armadillo_asynch/Test_24_model_2_async_overdamping_1.000000_prigid_1.000000/output")
video.add_frame("armadillo_ept5.png",caption_short_length)
if not testing_titles: video.add_directory("../render/armadillo_asynch/Test_24_model_2_async_overdamping_1.000000_prigid_0.500000/output")

video.add_frame("sphere1_bad_explicit.png",caption_long_length)
if not testing_titles: video.add_directory("../render/sphere_1_bad_explicit/Test_24_model_1_explicit_overdamping_1.000000_spheres_1_makeitbad_0.001000/output")
#video.add_frame("sphere1_bad_implicit.png",caption_length)
#if not testing_titles: video.add_directory("../render/sphere_1_bad_implicit/Test_24_model_1_implicit_overdamping_1.000000_spheres_1_makeitbad_0.001000/output",0,134)
video.add_frame("sphere1_bad_async_rigidpt5.png",caption_long_length)
if not testing_titles: video.add_directory("../render/sphere_1_bad_asynch/Test_24_model_1_async_overdamping_1.000000_spheres_1_prigid_0.500000_makeitbad_0.001000/output",0,120)
video.add_frame("title.png",title_length)

#if not testing_titles: video.add_directory("/output")

video.make_movie('asynchronous')

