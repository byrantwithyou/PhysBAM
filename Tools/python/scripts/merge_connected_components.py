#!/usr/bin/python

import os
import sys
import optparse
from struct import *

from physbam import *

#body="/home/mlentine/Desktop/mesh_test_6548.tri.gz"
#head="/home/mlentine/Desktop/mesh_test_9167.tri.gz"

parser = optparse.OptionParser("usage: %prog [options] input.tri output.tri")
parser.add_option('-d','--delete_connected',default=-1,type='int',help='delete connected components')
(options,args)=parser.parse_args()
if len(args)!=2: parser.error("invalid number of arguments")

tri=TRIANGULATED_SURFACE_f.Create()
Read_From_File("float",args[0],tri)
for i in range(1,len(args)):
    surface=TRIANGULATED_SURFACE_f.Create()
    Read_From_File("float",args[i],surface)
    for i in range(1,len(surface.mesh.elements)+1):
        tri.mesh.elements.Append(surface.mesh.elements[i])

def elem_compare(x,y):
    if(x[1][1]==y[1][1]):
        if(x[1][2]==y[1][2]): return x[1][3]-y[1][3]
        else: return x[1][2]-y[1][2]
    else: return x[1][1]-y[1][1]

def distance(x,y):
    return (x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])+(x[3]-y[3])*(x[3]-y[3]); 

def merge(x,y):
    for i in range(len(tri.mesh.elements),0,-1):
        if(tri.mesh.elements[i][1]==x): tri.mesh.elements[i][1]=y;
        if(tri.mesh.elements[i][2]==x): tri.mesh.elements[i][2]=y;
        if(tri.mesh.elements[i][3]==x): tri.mesh.elements[i][3]=y;
        if(tri.mesh.elements[i][1]==tri.mesh.elements[i][2] or tri.mesh.elements[i][1]==tri.mesh.elements[i][3] or tri.mesh.elements[i][2]==tri.mesh.elements[i][3]):
            tri.mesh.elements.Remove_Index_Lazy(i)

elements=[(i,tri.mesh.elements[i]) for i in range(1,len(tri.mesh.elements)+1)];
elements.sort(elem_compare);
for i in range(len(elements)-1):
    if(elements[i][1][1]==elements[i+1][1][1] and elements[i][1][2]==elements[i+1][1][2] and elements[i][1][3]==elements[i+1][1][3]):
        tri.mesh.elements.Remove_Index_Lazy(elements[i][0])

particle_list=[]
status=-1
for i in range(1,len(tri.particles)+1):
    particle_list.append((i,tri.particles.X[i]));
for i in range(len(tri.particles)-1):
    min_distance=distance(particle_list[i][1],particle_list[i+1][1]);
    min_index=i+1;
    for j in range(i+2,len(tri.particles)):
        new_distance=distance(particle_list[i][1],particle_list[j][1]);
        if(new_distance<min_distance):
            min_distance=new_distance;
            min_index=j;
    tmp=particle_list[i+1]
    particle_list[i+1]=particle_list[min_index]
    particle_list[min_index]=tmp

    new_status=(int(float(i)/float(len(tri.particles))*100))*1
    if (new_status is not status): print "%i"%new_status
    status=new_status

threshold=1e-8
for i in range(len(tri.particles)-1):
    print i
    if(distance(particle_list[i][1],particle_list[i+1][1])<threshold):
        merge(particle_list[i][0],particle_list[i+1][0])

Write_To_File("float","merged.tri.gz",tri);
