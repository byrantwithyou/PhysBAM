#!/usr/bin/python

import os
import sys
import optparse
from struct import *

from physbam import *

parser = optparse.OptionParser("usage: %prog [options] input.tri output.tri")
parser.add_option('-x','--x',default=1000,type='float',help='x plane to cut')
parser.add_option('-y','--y',default=1000,type='float',help='y plane to cut')
parser.add_option('-z','--z',default=1000,type='float',help='z plane to cut')
parser.add_option('-t','--tri',action="store_true",default=False,help='use a triangulated surface')
(options,args)=parser.parse_args()
if len(args)!=2: parser.error("invalid number of arguments")
 
tri=None
if(options.tri):tri=TRIANGULATED_SURFACE_f.Create()
else: tri=TETRAHEDRALIZED_VOLUME_f.Create()
Read_From_File("float",args[0],tri)

for i in range(len(tri.mesh.elements),0,-1):
    p1=tri.mesh.elements[i][1]
    p2=tri.mesh.elements[i][2]
    p3=tri.mesh.elements[i][3]
    if(tri.particles.X[p1][1]>options.x or tri.particles.X[p1][2]>options.y or tri.particles.X[p1][3]>options.z
        or tri.particles.X[p2][1]>options.x or tri.particles.X[p2][2]>options.y or tri.particles.X[p2][3]>options.z
        or tri.particles.X[p3][1]>options.x or tri.particles.X[p3][2]>options.y or tri.particles.X[p3][3]>options.z):
        tri.mesh.elements.Remove_Index_Lazy(i)

array=A_i();
tri.Discard_Valence_Zero_Particles_And_Renumber(array);
Write_To_File("float",args[1],tri)
