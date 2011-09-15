######################################################################
# Copyright 2007-2008, Geoffrey Irving, Andrew Selle.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# Python part of physbam python bindings
######################################################################
import sys
try:
    import dl
except ImportError:
    import DLFCN as dl
_dlopenflags = sys.getdlopenflags()
sys.setdlopenflags(dl.RTLD_NOW | dl.RTLD_GLOBAL)

import libphysbam_python_internal as _physbam
from libphysbam_python_internal import *

# restore dlopen flags
sys.setdlopenflags(_dlopenflags)

######################################################################
# Read/Write
######################################################################

class float(float):
    name='float'

class double(float):
    name='double'

def Get_Stream_Type(stream_type_name):
    try:
        stream_type_name=stream_type_name.name
    except:
        pass
    if stream_type_name=="float": return stream_type_float
    elif stream_type_name=="double": return stream_type_double
    else: raise TypeError("Invalid stream type '%s', should be 'float' or 'double'"%stream_type_name)

def Write_To_File(stream_type_name,filename,*objects):
    stream_type=Get_Stream_Type(stream_type_name)
    output=TYPED_OSTREAM(Safe_Open_Output(filename,True,True),stream_type)
    for i in objects:
        Write_Binary(output,i)

def Read_From_File(stream_type_name,filename,*objects):
    stream_type=Get_Stream_Type(stream_type_name)
    input=TYPED_ISTREAM(Safe_Open_Input(filename,True),stream_type)
    for i in objects:
        Read_Binary(input,i)

######################################################################
# Template constructors
######################################################################

def VECTOR(T,d): # T=float
    if isinstance(T,int): T,d=float,T
    return _physbam.__dict__['V%s%d'%(T.__name__[0],d)]

def MATRIX(T,m=None,n=None): # T=float
    if isinstance(T,int): T,m,n=float,T,m
    if n is None: n=m
    if m==n: return _physbam.__dict__['MATRIX_%s%d'%(T.__name__[0],m)]
    else: return _physbam.__dict__['MATRIX_%s%d%d'%(T.__name__[0],m,n)]

def DIAGONAL_MATRIX(T,m=None): # T=float
    if isinstance(T,int): T,m=float,T
    if m==1: return MATRIX(T,1)
    return _physbam.__dict__['DIAGONAL_MATRIX_%s%d'%(T.__name__[0],m)]

def SYMMETRIC_MATRIX(T,m=None): # T=float
    if isinstance(T,int): T,m=float,T
    if m==1: return MATRIX(T,1)
    return _physbam.__dict__['SYMMETRIC_MATRIX_%s%d'%(T.__name__[0],m)]

def UPPER_TRIANGULAR_MATRIX(T,m=None): # T=float
    if isinstance(T,int): T,m=float,T
    if m==1: return MATRIX(T,1)
    return _physbam.__dict__['UPPER_TRIANGULAR_MATRIX_%s%d'%(T.__name__[0],m)]

def SPHERE(center,radius):
    return _physbam.__dict__['SPHERE_'+type(center).__name__](center,radius)

def SOBOL(box):
    return _physbam.__dict__['SOBOL_'+type(box).__name__[4:]](box)

def BOX(TV):
    return _physbam.__dict__['BOX_'+TV.__name__]

def GRID(TV):
    return _physbam.__dict__['GRID_'+TV.__name__]

######################################################################
# Box iterators
######################################################################

class BOX_ITERATOR(object):
    def __init__(self,box):
        TV=type(box.min_corner)
        self.start=box.min_corner
        self.stop=box.max_corner
        self.d=len(self.start)
        self.current=TV(self.start)
        self.current[self.d]-=1 # adjust so that next() returns min_corner

    def next(self):
        d=self.d
        self.current[d]+=1
        if self.current[d]>self.stop[d]:
            for i in range(d,1,-1):
                self.current[i]=self.start[i]
                self.current[i-1]+=1
                if self.current[i-1]<=self.stop[i-1]:
                    break
            else:
                raise StopIteration
        return self.current

for d in 1,2,3:
    BOX(VECTOR(int,d)).__iter__=lambda box:BOX_ITERATOR(box)

######################################################################
# Miscellaneous functions
######################################################################

def sqr(x):
    return x*x

def cube(x):
    return x*x*x
