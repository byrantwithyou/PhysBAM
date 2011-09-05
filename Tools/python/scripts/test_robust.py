#!/usr/bin/python
#from physbam import *
from physbam import *

def stroke(psfile,simplex,fill):
    psfile.write("newpath\n")
    psfile.write("%f %f moveto\n"%simplex[0])
    for i in range(1,len(simplex)):
        psfile.write("%f %f lineto\n"%simplex[i])
    psfile.write(".01 setlinewidth\n")
    psfile.write("closepath\n")
    if fill:
        psfile.write(" .8 setgray\n")
        psfile.write("fill 0 setgray \n")
    else:
        psfile.write("stroke\n")

def draw(psfile,translate,tri,seg,status):
    psfile.write("gsave\n")
    psfile.write("%f %f translate\n"%translate)
    stroke(psfile,tri,True)
    stroke(psfile,seg,False)
    psfile.write(".01 .01  scale\n")
    psfile.write("0 0 moveto (%s) show\n"%repr(status))
    psfile.write("grestore\n")

    
def test(psfile,triangle_python,segment_python,i):
    triangle=V_Vf2_3(*map(lambda x: Vf2(*x),triangle_python))
    segment=V_Vf2_2(*map(lambda x: Vf2(*x),segment_python))
    
    #triangle=V_Vf2_3(Vf2(0,0),Vf2(0,1),Vf2(1,0))
    #segment=V_Vf2_2(Vf2(.25,.25),Vf2(1.25,0))
    interactions=ROBUST_SIMPLEX_INTERACTIONS_Vf2()
    robust=False
    intersecting,robust=interactions.Intersection_Test(triangle,segment)


    draw(psfile,(0,i),triangle_python,segment_python,(intersecting,robust))
    
    print "intersecting %d, robust %d"%(intersecting,robust)
    
psfile=open("test_robust.ps","w")
psfile.write("%!PS\n")
psfile.write("50 50 translate\n")
psfile.write("""/Helvetica findfont
10 scalefont
setfont
""")
    
psfile.write("100 100 scale\n")

triangle_python=((0.,0.),(0.,1.),(1.,0.))
segment_python=((2.25,.25),(4.25,0.))
test(psfile,triangle_python,segment_python,0)

triangle_python=((0.,0.),(0.,1.),(1.,0.))
segment_python=((0.25,.25),(2.25,0.))
test(psfile,triangle_python,segment_python,1)

triangle_python=((0.,0.),(0.,1.),(1.,0.))
segment_python=((0.25,.25),(2.25,0.25))
test(psfile,triangle_python,segment_python,2)

triangle_python=((0.,0.),(0.,1.),(1.,0.))
segment_python=((1.25,.25),(2.25,0.25))
test(psfile,triangle_python,segment_python,3)

triangle_python=((1.,0.),(0.,1.),(1.,0.))
segment_python=((0,1),(1,0))
test(psfile,triangle_python,segment_python,4)


test(psfile,((0.25,.65),(5,.5),(-.1,0.)),((2.25,0.25),(-.5,.1)),5)
test(psfile,((0.25,.65),(5,.5),(-.1,0.)),((2.25,0.25),(2.5,.7)),6)

psfile.write("showpage\n")
