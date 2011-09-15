#!/usr/bin/python
######################################################################
# Copyright 2008, Geoffrey Irving.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# Test shape.py
######################################################################

import unittest
import tempfile
import shutil
import copy

from physbam import *
from math import *

T=float
TV=Vf3

class shape_tests(unittest.TestCase):

    def setUp(self):
        self.shapes=[
            SPHERE_Vf3(TV(),1),
            BOX_Vf3((-1,-2,-3),(1,2,3)),
            PLANE_f((1,0,0),(0,1,0),(0,0,1)),
            CYLINDER_f((0,-.5,0),(0,.5,0),1),
            RING_f((0,-.5,0),(0,.5,0),3,2)]

    def test_consistency(self): 
        n=1000
        for shape in self.shapes:
            print shape
            if isinstance(shape,BOX_Vf3):
                box=shape
            elif isinstance(shape,PLANE_f):
                box=BOX_Vf3.Unit_Box()
            else:
                box=shape.Bounding_Box()
            scale=max(box.Edge_Lengths())
            sobol=SOBOL(box.Thickened(scale/pi))
            inner_box=BOX_Vf3.Empty_Box()
            for _ in range(n):
                X=sobol.Get_Vector()
                if X in shape:
                    inner_box.Enlarge_To_Include_Point(X)
                phi=shape.Signed_Distance(X)
                normal=shape.Normal(X)
                project=X-phi*normal
                phi2=shape.Signed_Distance(project)
                assert abs(phi2)<scale/100,'%s = (%s) - %s * (%s), %s'%(project,X,phi,normal,phi2)
                assert shape.Inside(X,0)==shape.Lazy_Inside(X)==(phi<=0) # TODO: nonzero thickness
                assert shape.Outside(X,0)==shape.Lazy_Outside(X)==(phi>0) # TODO: nonzero thickness
                surface=shape.Surface(X)
                assert abs(shape.Signed_Distance(surface))<scale/100
                surface_phi=(surface-X).Magnitude()
                assert abs(surface_phi-abs(phi))<scale/100,'%s != %s'%(surface_phi,phi)
                # TODO: test Boundary and Principal_Curvatures
            if not isinstance(shape,PLANE_f):
                assert inner_box.min_corner in box
                assert inner_box.max_corner in box
                assert max(box.Edge_Lengths()-inner_box.Edge_Lengths())<scale/4,'%s != %s'%(box,inner_box)

    def test_generate_triangles(self):
        tolerance=1e-5
        for shape in self.shapes:
            print shape
            #surface=shape.Generate_Triangles()
            if not isinstance(shape,PLANE_f):
                assert not surface.mesh.Non_Manifold_Nodes()
            assert surface.mesh.Orientations_Consistent()
            particles=surface.particles
            for X in particles.X:
                assert abs(shape.Signed_Distance(X))<tolerance
            for t in range(1,len(surface.mesh.elements)+1):
                i,j,k=surface.mesh.elements[t]
                X=(particles.X[i]+particles.X[j]+particles.X[k])/3
                shape_normal=shape.Normal(X)
                surface_normal=surface.Normal(X,t)
                dot=TV.Dot_Product(shape_normal,surface_normal)
                assert dot>.7,'(%s) . (%s) = %s'%(shape_normal,surface_normal,dot)
                
######################################################################
# Main
######################################################################

if __name__ == '__main__':
    unittest.main()
