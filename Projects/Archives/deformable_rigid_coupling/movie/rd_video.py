#!/usr/bin/python

from video import *
import sys

header_color='#ffaa00'
subtle_color='#808080'

slide("title.png",40,[(40,"Two-way Coupling"),(40,"of Rigid"),("and Deformable Bodies"),(40,""),(36,"SCA 2008 #1013")])

slide("torsion-spring.png",24,[(36,"Deformable bar",header_color),"","embeddings, joint constraints,","contact, and collisions","deformable bar (blue)","rigid boxes (silver)"])
slide("row-spheres-far.png",24,[(36,"Impact propagation",header_color),"","rigid balls (silver)","incompressible deformable balls (yellow)"])
slide("friction.png",24,[(36,"Friction test",header_color),"","analytic solution (white)","rigid box (silver)","deformable box (orange)","particle (blue)"])
slide("cube-stack.png",24,[(36,"Stack",header_color),"","rigid bodies (silver)","deformable boxes (orange)"])
slide("trampoline.png",24,[(36,"Trampoline",header_color),"","cloth (white)","rigid balls (black)","deformable tori (colored)"])

slide("stress.png",32,["Deformable and rigid bodies","with contact and constraints","are robustly simulated","even under extreme scenarios."])
slide("twisting-chain.png",24,[(36,"Chain under high stress",header_color),"","rigid tori (silver)","articulated rigid loops (rainbow)","deformable tori (yellow)"])
slide("scaling.png",36,["Our method scales","to large numbers","of bodies."])
slide("ring-drop.png",24,[(36,"Large pile",header_color),"","320 deformable objects","1280 rigid bodies","(1600 objects)"])

# oscopy.
slide("creatures1.png",36,["Two-way coupling is","important for modeling","life-like creatures","that interact","with their environment."])
slide("creatures2.png",36,["Each creature","is modeled with","an embedded PD-controlled","articulated rigid skeleton","and a deformable exterior."])
slide("snake.png",24,[(36,"Sidewinding snake",header_color),"","12 joints","13 rigid bodies","deformable exterior","(semitransparent)"])
slide("snake-view2.png",36,["Snake climbs stairs"])
slide("snake-view1.png",36,["Another view","with snake","rendered opaque"])

slide("maggot.png",24,[(36,"Wriggling Maggot",header_color),"","2 joints","3 rigid bodies","deformable exterior"])
slide("maggot-trapped.png",36,["Maggot trapped","under a rigid torus"])
slide("maggot-bowl.png",36,["20 maggots","dropped","into a bowl"])
slide("maggot-rings.png",36,["... with 20 rigid rings"])

slide("fish.png",24,[(36,"Flopping fish",header_color),"","4 joints","5 rigid bodies","deformable exterior"])

slide("end.png",40,[(40,"Two-way Coupling"),(40,"of Rigid"),("and Deformable Bodies"),(40,""),(36,"SCA 2008 #1013")])

slide("closer.png",36,["Closer"])
slide("again.png",36,["Again"])
slide("slower.png",36,["Slower"])

video_dir="rd_video_tmp-100"
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

video.add_frame("stress.png",caption_longer_length)
video.add_frame("twisting-chain.png",caption_longer_length)
if not testing_titles: video.add_directory("../render/twisting-chain-reverse/output-final")

video.add_frame("cube-stack.png",caption_length)
if not testing_titles: video.add_directory("../render/stack6/output-final")

video.add_frame("row-spheres-far.png",caption_long_length)
if not testing_titles: video.add_directory("../render/row-sphere-incomp-hires-far/output-final")
video.add_frame("closer.png",caption_tinier_length)
if not testing_titles: video.add_directory("../render/row-sphere-incomp-hires-near/output-final")
video.add_frame("slower.png",caption_tinier_length)
if not testing_titles: video.add_directory("row-sphere-incomp-hires-near-slow")

video.add_frame("torsion-spring.png",caption_longer_length)
if not testing_titles: video.add_directory("torsion-spring")

video.add_frame("trampoline.png",caption_long_length)
if not testing_titles: video.add_directory("../render/trampoline4/output-final")

video.add_frame("scaling.png",caption_length)
video.add_frame("ring-drop.png",caption_long_length)
if not testing_titles: video.add_directory("ring-drop")

video.add_frame("friction.png",caption_longer_length)
if not testing_titles: video.add_directory("cubes-friction")
if not testing_titles: video.add_frame("cubes-friction/cubes-friction.00096.png",caption_tiny_length)

video.add_frame("creatures1.png",caption_long_length)
video.add_frame("creatures2.png",caption_long_length)
video.add_frame("snake.png",caption_long_length)
if not testing_titles: video.add_directory("../render/snake-course-close/output-final")
video.add_frame("snake-view2.png",caption_short_length)
if not testing_titles: video.add_directory("../render/snake-course-view2/output-final")
# video.add_frame("snake-view1.png",caption_long_length)
# if not testing_titles: video.add_directory("../render/snake-course-view1/output-final")

video.add_frame("maggot.png",caption_length)
if not testing_titles: video.add_directory("../render/maggot2/output-final")
video.add_frame("maggot-trapped.png",caption_length)
if not testing_titles: video.add_directory("trapped") # no preroll
video.add_frame("maggot-bowl.png",caption_length)
if not testing_titles: video.add_directory("../render/bowl-20/output-final")
video.add_frame("maggot-rings.png",caption_length)
if not testing_titles: video.add_directory("../render/maggot-rings/output-final")

video.add_frame("fish.png",caption_length)
if not testing_titles: video.add_directory("../render/floppy-fish/output-final")

video.add_frame("end.png",title_length)

video.make_movie('rigid-deformable')

