#!/usr/bin/python

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

def draw(psfile,translate,seg1,seg2,status):
    psfile.write("gsave\n")
    psfile.write("%f %f translate\n"%translate)
    stroke(psfile,seg1,False)
    stroke(psfile,seg2,False)
    psfile.write(".005 .005  scale\n")
    psfile.write("0 0 moveto (%s) show\n"%status)
    psfile.write("grestore\n")

def Test(psfile,segment1_python,segment2_python,count):
    segment1=V_Vf2_2(*map(lambda x: Vf2(*x),segment1_python))
    segment2=V_Vf2_2(*map(lambda x: Vf2(*x),segment2_python))

    s=SIMPLEX_INTERACTIONS_f()
    #SIMPLEX_INTERACTIONS_f.Test()
    #s.Test()
    w1,w2=Vf1(),Vf1()
    intersected=s.Two_Segment_Intersection_Barycentric_Coordinates(segment1,segment2,w1,w2)
    draw(psfile,(count/10*2,count%10),segment1_python,segment2_python,"intersected=%d w1=%f w2=%f"%(intersected,w1.x,w2.x))
    
psfile=open("test_robust.ps","w")
psfile.write("%!PS\n")
psfile.write("50 50 translate\n")
psfile.write("""/Helvetica findfont
10 scalefont
setfont
""")
    
    
psfile=open("test_robust.ps","w")
psfile.write("%!PS\n")
psfile.write("50 50 translate\n")
psfile.write("""/Helvetica findfont
10 scalefont
setfont

50 50 scale
""")

counter=0
def count():
    global counter
    counter+=1
    return counter

Test(psfile,((0,0),(1,0)),((.5,-.5),(.5,.5)),count())
Test(psfile,((0,0),(1,0)),((0,0),(.5,.5)),count())
Test(psfile,((0,0),(1,0)),((0,1),(1,1)),count())
import random
r=random.Random()
for i in range(0,20):
    def u():
        return r.uniform(0,1)
    Test(psfile,((u(),u()),(u(),u())),((u(),u()),(u(),u())),count())

psfile.write("showpage\n")


