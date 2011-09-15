from video import *
import sys

slide("title.png",40,[(36,"paper_0432"),(40,""),(40,"Volume Conserving"),(40,"Finite Element Simulations"),(40,"of Deformable Models")])

slide("sphere-drop-poisson.png",36,["Increasing Poisson's ratio","increases incompressibility","but leads to","nonphysical locking"])
slide("sphere-drop-3.png",36,["Poisson's ratio .3","","40% max volume loss"])
slide("sphere-drop-4.png",36,["Poisson's ratio .4","","26% max volume loss"])
slide("sphere-drop-45.png",36,["Poisson's ratio .45","","16% max volume loss"])
slide("sphere-drop-49.png",36,["Poisson's ratio .49","","4.9% max volume loss","locking"])
slide("sphere-drop-499.png",36,["Poisson's ratio .499","","2.2% max volume loss","severe locking"])
slide("sphere-drop-5.png",36,["Our method","","<1% max volume loss","no locking"])

slide("sphere-drop-stiffness.png",36,["Our method","","varying stiffness"])
slide("sphere-drop-low.png",36,["Low stiffness","<1% max volume loss"])#note volume loss - skip a line first, for all three
slide("sphere-drop-medium.png",36,["Medium stiffness","<1% max volume loss"])
slide("sphere-drop-high.png",36,["High stiffness","<1% max volume loss"])

slide("sphere-stairs.png",36,["Incompressible sphere","falls down stairs"])

slide("armadillo-drop-title.png",36,["Incompressible armadillo","drops onto ground"])
slide("armadillo-drop-low.png",36,["Low stiffness","<1% max volume loss","(112k tetrahedra)"])
slide("armadillo-drop-hi.png",36,["Higher stiffness","<1% max volume loss","(112k tetrahedra)"])

slide("armadillo-stairs-title.png",36,["Incompressible armadillo","falls down stairs"])
slide("armadillo-stairs-low.png",36,["Low stiffness","<1% max volume loss","(112k tetrahedra)"])
slide("armadillo-stairs-medium.png",36,["Higher stiffness","<1% max volume loss","(112k tetrahedra)"])
slide("armadillo-stairs-high.png",36,["Even higher stiffness","<1% max volume loss","(112k tetrahedra)"])

slide("tori-pile-title-100.png",36,["100 incompressible tori","falling into a pile","(1.25m tetrahedra)"])

slide("end.png",40,[(36,"paper_0432"),(40,""),(40,"Volume Conserving"),(40,"Finite Element Simulations"),(40,"of Deformable Models")])

slide("again.png",36,["Again"])
slide("slower.png",36,["Slower"])

video_dir="cloth_video_tmp-100"
if os.path.isdir(video_dir):
    print "%s already exists... delete? (y/n) "%video_dir,
    c=sys.stdin.readline()
    if c=="y\n": shutil.rmtree(video_dir)
    else: sys.exit(1)

video=VIDEO(video_dir)

testing_titles=False
if testing_titles:
    title_length=video.fps*5/10
    caption_length=int(video.fps*3.5/10)
    caption_long_length=int(video.fps*5/10)
    caption_short_length=int(video.fps*2/10)
    caption_tiny_length=int(video.fps*1/10)
else:
    title_length=video.fps*5*1.2
    caption_length=int(video.fps*3.5*1.2)
    caption_long_length=int(video.fps*5*1.2)
    caption_short_length=int(video.fps*2*1.2)
    caption_tiny_length=int(video.fps*1*1.2)
#testing_titles=True
print caption_length


video.add_frame("title.png",title_length)

video.add_frame("sphere-drop-poisson.png",caption_long_length)
video.add_frame("sphere-drop-3.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-300")
video.add_frame("sphere-drop-4.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-400")
video.add_frame("sphere-drop-45.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-450")
video.add_frame("sphere-drop-49.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-490")
video.add_frame("sphere-drop-499.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-499")
video.add_frame("sphere-drop-5.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-500")

video.add_frame("sphere-drop-stiffness.png",caption_long_length)
video.add_frame("sphere-drop-low.png",caption_short_length)
if not testing_titles:
    video.add_directory("stiff_0.5")
video.add_frame("sphere-drop-medium.png",caption_short_length)
if not testing_titles:
    video.add_directory("stiff_3")
video.add_frame("sphere-drop-high.png",caption_short_length)
if not testing_titles:
    video.add_directory("stiff_20")

video.add_frame("sphere-stairs.png",caption_length)
if not testing_titles:
    video.add_directory("stair-good")

video.add_frame("armadillo-drop-title.png",caption_length)
video.add_frame("armadillo-drop-low.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo2")
video.add_frame("slower.png",caption_tiny_length)
if not testing_titles:
    video.add_directory("armadillo2-slow")
video.add_frame("armadillo-drop-hi.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo5")
video.add_frame("slower.png",caption_tiny_length)
if not testing_titles:
    video.add_directory("armadillo5-slow")

video.add_frame("armadillo-stairs-title.png",caption_length)
video.add_frame("armadillo-stairs-low.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-2")
video.add_frame("armadillo-stairs-medium.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-4")
video.add_frame("armadillo-stairs-high.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-8")

video.add_frame("tori-pile-title-100.png",caption_length)
if not testing_titles:
    video.add_directory("tori-100")
video.add_frame("slower.png",caption_tiny_length)
if not testing_titles:
    video.add_directory("tori-100-dup")

video.add_frame("end.png",title_length)


video.make_movie('incompressible')





video=VIDEO(video_dir+'p')
video.add_frame("sphere-drop-poisson.png",caption_long_length)
video.add_frame("sphere-drop-3.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-300")
video.add_frame("sphere-drop-4.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-400")
video.add_frame("sphere-drop-45.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-450")
video.add_frame("sphere-drop-49.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-490")
video.add_frame("sphere-drop-499.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-499")
video.add_frame("sphere-drop-5.png",caption_length)
if not testing_titles:
    video.add_directory("poisson-500")
video.make_movie('incompressible-p')

video=VIDEO(video_dir+'s')
video.add_frame("sphere-drop-stiffness.png",caption_long_length)
video.add_frame("sphere-drop-low.png",caption_short_length)
if not testing_titles:
    video.add_directory("stiff_0.5")
video.add_frame("sphere-drop-medium.png",caption_short_length)
if not testing_titles:
    video.add_directory("stiff_3")
video.add_frame("sphere-drop-high.png",caption_short_length)
if not testing_titles:
    video.add_directory("stiff_20")
video.make_movie('incompressible-s')

video=VIDEO(video_dir+'ss')
video.add_frame("sphere-stairs.png",caption_length)
if not testing_titles:
    video.add_directory("stair-good")
video.make_movie('incompressible-ss')

video=VIDEO(video_dir+'ad')
video.add_frame("armadillo-drop-title.png",caption_length)
video.add_frame("armadillo-drop-low.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo2")
video.add_frame("armadillo-drop-hi.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo5")
video.make_movie('incompressible-ad')

video=VIDEO(video_dir+'as')
video.add_frame("armadillo-stairs-title.png",caption_length)
video.add_frame("armadillo-stairs-low.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-2")
video.add_frame("armadillo-stairs-medium.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-4")
video.add_frame("armadillo-stairs-high.png",caption_length)
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-8")
video.make_movie('incompressible-as')

video=VIDEO(video_dir+'t')
video.add_frame("tori-pile-title-100.png",caption_length)
if not testing_titles:
    video.add_directory("tori-100")
video.make_movie('incompressible-t')



video=VIDEO(video_dir+'p-3')
if not testing_titles:
    video.add_directory("poisson-300")
video.make_movie('incompressible-p-3')
video=VIDEO(video_dir+'p-4')
if not testing_titles:
    video.add_directory("poisson-400")
video.make_movie('incompressible-p-4')
video=VIDEO(video_dir+'p-45')
if not testing_titles:
    video.add_directory("poisson-450")
video.make_movie('incompressible-p-45')
video=VIDEO(video_dir+'p-49')
if not testing_titles:
    video.add_directory("poisson-490")
video.make_movie('incompressible-p-49')
video=VIDEO(video_dir+'p-499')
if not testing_titles:
    video.add_directory("poisson-499")
video.make_movie('incompressible-p-499')
video=VIDEO(video_dir+'p-5')
if not testing_titles:
    video.add_directory("poisson-500")
video.make_movie('incompressible-p-5')

video=VIDEO(video_dir+'s-ll')
if not testing_titles:
    video.add_directory("stiff_0.5")
video.make_movie('incompressible-s-l')
video=VIDEO(video_dir+'s-m')
if not testing_titles:
    video.add_directory("stiff_3")
video.make_movie('incompressible-s-m')
video=VIDEO(video_dir+'s-h')
if not testing_titles:
    video.add_directory("stiff_20")
video.make_movie('incompressible-s-h')

video=VIDEO(video_dir+'ss-x')
if not testing_titles:
    video.add_directory("stair-good")
video.make_movie('incompressible-ss-x')

video=VIDEO(video_dir+'ad-l')
if not testing_titles:
    video.add_directory("armadillo2")
video.make_movie('incompressible-ad-l')
video=VIDEO(video_dir+'ad-h')
if not testing_titles:
    video.add_directory("armadillo5")
video.make_movie('incompressible-ad-h')

video=VIDEO(video_dir+'as-l')
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-2")
video.make_movie('incompressible-ss-l')
video=VIDEO(video_dir+'as-m')
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-4")
video.make_movie('incompressible-ss-m')
video=VIDEO(video_dir+'as-h')
if not testing_titles:
    video.add_directory("armadillo-stairs-stiffness-8")
video.make_movie('incompressible-ss-h')

video=VIDEO(video_dir+'t-x')
if not testing_titles:
    video.add_directory("tori-100")
video.make_movie('incompressible-t-x')

