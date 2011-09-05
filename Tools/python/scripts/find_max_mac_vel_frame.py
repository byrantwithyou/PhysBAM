#!/usr/bin/python
import sys
import physbam
import os
import math

u=physbam.ARRAYS_3D_f()
v=physbam.ARRAYS_3D_f()
w=physbam.ARRAYS_3D_f()

file=open("last_frame",'r');
if file:
    last_frame=int(file.readline());
file.close()

old_max=1;
for i in range(125,last_frame):
    mac_filename="mac_velocities.%d"%i
    physbam.Read_From_File("float",mac_filename,u,v,w)
    max=0.0
    m=u.m
    n=u.n
    mn=u.mn
    if(v.m<m): m=v.m
    if(w.m<m): m=w.m
    if(v.n<n): n=v.n
    if(w.n<n): n=w.n
    if(v.mn<mn): mn=v.mn
    if(w.mn<mn): mn=w.mn
    mn-=3
    n-=3
    m-=3
    for j in range(1,m+1):
        for k in range(1,n+1):
            for jk in range(1,mn+1):
                index=physbam.Vi3(j,k,jk);
                magnitude=math.sqrt(u.__getitem__(index)*u.__getitem__(index)+v.__getitem__(index)*v.__getitem__(index)+w.__getitem__(index)*w.__getitem__(index))
                if magnitude>max: max=magnitude
    print "Max at frame %d is %f"%(i,max)
    print "Jump is %f"%(max/old_max);
    old_max=max;
    if old_max==0: old_max=1;
