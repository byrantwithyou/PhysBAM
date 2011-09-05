#!/usr/bin/python
######################################################################
# Copyright 2007, Geoffrey Irving.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# Script sphere.py
######################################################################

from physbam import *
from math import *
import random
import collections
import unittest

epsilon=numeric_limits_f.epsilon()

def Simple_Bounding_Sphere(points):
    box=BOX(type(points[0])).Bounding_Box(points)
    center=box.Center()
    radius=sqrt(max((p-center).Magnitude_Squared() for p in points))
    sphere=SPHERE(center,radius)
    sphere.bad=True
    return sphere

def Assert_Inside(s1,s2):
    tolerance=2*epsilon*max(s1.center.Magnitude(),s1.radius)
    phi=s2.Signed_Distance(s1.center)+s1.radius
    try:
        assert(phi<=tolerance)
    except AssertionError:
        print 's1 %r, s2 %r, phi %r, tolerance %r'%(s1,s2,phi,tolerance)
        raise

def Check_Inclusion(sphere,points):
    large=Simple_Bounding_Sphere(points)
    tolerance=2*epsilon*max(sphere.center.Magnitude(),sphere.radius)
    max_phi=max(sphere.Signed_Distance(p) for p in points)
    large_phi=large.Signed_Distance(sphere.center)+sphere.radius
    try:
        assert(max_phi<=tolerance)
        assert(sphere.radius<=1.04*large.radius+tolerance)
    except AssertionError:
        print 'points:',points
        print 'sphere:',sphere
        print 'max_phi:',max_phi
        print 'large:',large
        print 'large_phi:',large_phi
        print 'tolerance %g'%tolerance
        raise

def Check_Sphere(sphere,support,points):
    Check_Inclusion(sphere,points)
    if 'bad' in sphere.__dict__: return # no support if we switched to approximate version
    assert(len(support)<=len(sphere.center)+1)
    support_error=max(abs(sphere.Signed_Distance(p)) for p in support)
    tolerance=10*pow(epsilon,1./len(sphere.center))*max(sphere.center.Magnitude(),sphere.radius)
    try:
        assert(support_error<=tolerance)
    except AssertionError:
        print 'points:',points
        print 'support:',support
        print 'sphere %r, support_error %g'%(sphere,support_error)
        print 'tolerance %g'%tolerance
        raise

def Approximate_Bounding_Sphere(points):
    # from Jack Ritter (ftp://ftp-graphics.stanford.edu/pub/Graphics/RTNews/html/rtnews7b.html#art4)
    box=BOX(type(points[0])).Bounding_Box(points)
    sphere=SPHERE(box.Center(),max(box.Edge_Lengths())/2)
    for p in points:
        DX=p-sphere.center
        sqr_distance=DX.Magnitude_Squared()
        if sqr_distance>sqr(sphere.radius):
            distance=sqrt(sqr_distance)
            shift=.5*(distance-sphere.radius)
            sphere=SPHERE(sphere.center+shift/distance*DX,sphere.radius+shift)
    sphere.bad=True
    Check_Inclusion(sphere,points)
    return sphere

def Interpolating_Sphere(points):
    if not points:
        return ()
    base,rest=points[0],points[1:]
    if not rest:
        return SPHERE(base,0)
    elif len(rest)==1:
        other,=rest
        return SPHERE((other+base)/2,(other-base).Magnitude()/2)
    Q=MATRIX(len(base),len(rest))(*(p-base for p in rest))
    A=Q.Normal_Equations_Matrix()
    b=A.Diagonal_Part().To_Vector()/2
    c=A.Robust_Solve_Linear_System(b) # may return garbage, but avoids exceptions
    if [() for x in c if x<0] or c.Sum()>1: # detect garbage
        return Approximate_Bounding_Sphere(points)
    relative_center=Q*c
    sqr_radius=relative_center.Magnitude_Squared()
    center=relative_center+base
    # enlarge radius to account for numerical error if necessary
    for p in rest:
        sqr_radius=max(sqr_radius,(p-center).Magnitude_Squared())
    sphere=SPHERE(center,sqrt(sqr_radius))
    # check whether approximate sphere enclosing box is better (only happens due to numerical error)
    simple_sphere=Approximate_Bounding_Sphere(points)
    if sphere.radius>simple_sphere.radius:
        sphere=simple_sphere
    Check_Inclusion(sphere,points)
    return sphere

bounding_sphere_count=0
simple_bounding_sphere_incremental_fix_count=0

def Bounding_Sphere(points,rand=random.random):
    global bounding_sphere_count
    bounding_sphere_count+=1

    assert(len(points)>=1)
    d=len(points[0])

    def BS(interior,boundary):
        if len(boundary)==d+1:
            sphere=Interpolating_Sphere(boundary)
            flag=False
            old_radius=sphere.radius
            for p in interior:
                if p not in sphere:
                    sphere.radius=(p-sphere.center).Magnitude()
                    flag=True
            if flag:
                print old_radius/sphere.radius
                global simple_bounding_sphere_incremental_fix_count
                simple_bounding_sphere_incremental_fix_count+=1
                simple_sphere=Approximate_Bounding_Sphere(boundary+tuple(interior))
                if simple_sphere.radius<sphere.radius:
                    sphere=simple_sphere
                sphere.bad=True
            return sphere,boundary
        active=collections.deque()
        sphere,support=Interpolating_Sphere(boundary),boundary
        while interior:
            p=interior.popleft()
            if p not in sphere:
                sphere,support=BS(active,boundary+(p,))
                active.appendleft(p) # move to front
            else:
                active.append(p)
        interior.extend(active) # put everything back in interior
        return sphere,support

    interior=collections.deque(points)
    random.shuffle(interior,random=rand)
    sphere,support=BS(interior,())
    simple_sphere=Approximate_Bounding_Sphere(points)
    if sphere.radius>simple_sphere.radius:
        sphere=simple_sphere
    Check_Sphere(sphere,support,points)
    return sphere

class SPHERE_TESTS(unittest.TestCase):
    def setUp(self):
        self.rand=RANDOM_NUMBERS_f()
        self.rand.Set_Seed(92323)

    def test_approximate(self):
        for d in [2,3]:
            for _ in range(100):
                for c in range(1,10)+[1000]:
                    points=[self.rand.Get_Number()*self.rand.Get_Direction(d) for i in range(c)]
                    Approximate_Bounding_Sphere(points)

    def test_approximate_cospherical(self):
        for d in [2,3]:
            for _ in range(100):
                for c in range(1,10):
                    points=[self.rand.Get_Direction(d) for i in range(c)]
                    Approximate_Bounding_Sphere(points)

    def test_small(self):
        for d in [2,3]:
            for _ in range(100):
                for c in range(1,10)+[1000]:
                    points=[self.rand.Get_Number()*self.rand.Get_Direction(d) for i in range(c)]
                    Bounding_Sphere(points,rand=self.rand.Get_Number)

    def test_cospherical(self):
        for d in [2,3]:
            for _ in range(100):
                for c in range(1,10):
                    points=[self.rand.Get_Direction(d) for i in range(c)]
                    Bounding_Sphere(points,rand=self.rand.Get_Number)

    def test_coplanar(self):
        for d in [2,3]:
            for _ in range(100):
                for c in range(1,10):
                    base=self.rand.Get_Direction(d)
                    normal=self.rand.Get_Direction(d) 
                    points=[base+self.rand.Get_Direction(d).Projected_Orthogonal_To_Unit_Direction(normal) for i in range(c)]
                    Bounding_Sphere(points,rand=self.rand.Get_Number)

if __name__ == '__main__':
    try:
        unittest.main()
    finally:
        print 'bounding_sphere_count',bounding_sphere_count
        print 'simple_bounding_sphere_incremental_fix_count',simple_bounding_sphere_incremental_fix_count

