#!/usr/bin/python

import os
import sys
import optparse
from struct import *

from physbam import *

parser = optparse.OptionParser("usage: %prog [options] input.tri")
(options,args)=parser.parse_args()
if len(args)!=1: parser.error("invalid number of arguments")

tri=TRIANGULATED_SURFACE_f.Create()
Read_From_File("float",args[0],tri)

delete=[]
for i in range(1,len(tri.particles)+1):
    found=False
    for j in range(1,len(tri.mesh.elements)+1):
        p1=tri.mesh.elements[j][1]
        p2=tri.mesh.elements[j][1]
        p3=tri.mesh.elements[j][1]
        if(p1==i or p2==i or p3==i):
            found=True
    if(not found): delete.append(i)

print delete
