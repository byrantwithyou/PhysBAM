#!/usr/bin/python
######################################################################
# Copyright 2007, Geoffrey Irving.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# Script sphere_hierarchy.py
######################################################################

from physbam import *
from math import *
from sphere import Bounding_Sphere
import os
import random
import unittest

epsilon=numeric_limits_f.epsilon()

def Sobol_Integrate(f,box,tolerance,epoch=1000):
    sobol=SOBOL(box)
    n=0
    sum=0.
    while True:
        min_mean,max_mean=1e10,-1e10
        for _ in range(epoch):
            n+=1
            y=f(sobol.Get_Vector())
            sum+=y
            mean=sum/n
            min_mean=min(min_mean,mean)
            max_mean=max(max_mean,mean)
        if max_mean-min_mean<=.5*tolerance*abs(sum)/n:
            return box.Size()*sum/n

def Sphere_Intersection_Volume(s1,s2):
    if s1.radius>s2.radius: s1,s2=s2,s1 # smaller sphere first
    d=(s1.center-s2.center).Magnitude()
    r1,r2=s1.radius,s2.radius
    if d>=r1+r2: return 0 # no overlap
    elif d+r1<=r2: return s1.Size()
    else: return pi*sqr(r1+r2-d)*(sqr(d)+2*d*(r1+r2)-3*sqr(r1-r2))/(12*d)

class SPHERE_HIERARCHY(object):
    def __init__(self,particles):
        self.particles=particles
        self.TV=Vf3

    class NODE(object):
        def __init__(self,X=None):
            if X:
                self.X=X
                self.sphere=SPHERE(X,0)
            self.children=[]

        def Points(self):
            if self.children:
                return sum((c.Points() for c in self.children),())
            else:
                return self.X,

        def Slow_Update(self):
            self.sphere=Bounding_Sphere(self.Points())

    def Initialize_From_KD_Tree(self,particles_per_group):
        if particles_per_group==1: particles_per_group=0
        hierarchy=PARTICLE_HIERARCHY_Vf3(self.particles,True,particles_per_group)
        X=self.particles.X
        def Convert(box):
            node=self.NODE()
            if box<=hierarchy.leaves:
                if hierarchy.particles_per_group:
                    node.children.extend(self.NODE(X[p]) for p in hierarchy.particles_in_group[box])
                else:
                    node.X=X[box]
            else:
                node.children.extend(Convert(b) for b in hierarchy.children[box-hierarchy.leaves])
            return node
        self.root=Convert(hierarchy.root)

    def Incremental_Initialize(self):
        def Add(X):
            best_slot=[(1e10,None)]
            def Find(node,parent,base_cost):
                if not node.children: # we would need to create a new sphere containing the two points in question
                    sphere=Bounding_Sphere([node.X,X])
                    size=Sphere_Intersection_Volume(sphere,parent)
                    best_slot[0]=min(best_slot[0],(base_cost+2*size,node))
                else:
                    if best_slot[0][0]<=base_cost:
                        return
                    sphere=node.sphere
                    size=Sphere_Intersection_Volume(sphere,parent)
                    if X not in sphere:
                        old_size=size 
                        sphere=Bounding_Sphere(node.Points()+(X,))
                        size=Sphere_Intersection_Volume(sphere,parent)
                        # add cost due to expanding this sphere to include X
                        base_cost+=(size-old_size)*len(node.children)
                        for c in node.children:
                            old_child_size=Sphere_Intersection_Volume(c.sphere,node.sphere)
                            new_child_size=Sphere_Intersection_Volume(c.sphere,sphere)
                            base_cost+=(new_child_size-old_child_size)*len(node.children)
                    # consider adding X to this node
                    best_slot[0]=min(best_slot[0],(base_cost+size,node))
                    # consider adding X to each child
                    for c in node.children:
                        if X in c.sphere:
                            Find(c,sphere,base_cost)
                    for c in node.children:
                        if X not in c.sphere:
                            Find(c,sphere,base_cost)
            Find(self.root,self.root.sphere,0)
            target=best_slot[0][1]
            def Insert(node):
                if node==target:
                    if node.children:
                        node.children.append(self.NODE(X))
                        if X not in node.sphere:
                            node.sphere=Bounding_Sphere(node.Points())
                    else:
                        node.children.append(self.NODE(node.X))
                        node.children.append(self.NODE(X))
                        node.sphere=Bounding_Sphere((node.X,X))
                        del node.X
                    return True
                else:
                    for c in node.children:
                        if Insert(c):
                            if X not in node.sphere:
                                node.sphere=Bounding_Sphere(node.Points())
                            return True
                    return False
            Insert(self.root)
        points=list(self.particles.X)
        random.shuffle(points)
        X=points.pop()
        self.root=self.NODE(X)
        for X in points:
            Add(X)
        assert(self.Leaf_Count()==len(self.particles.X))

    def Slow_Update(self):
        def Update(node):
            if node.children:
                points=sum((Update(c) for c in node.children),())
                node.sphere=Bounding_Sphere(points)
                return points
            else:
                node.sphere=SPHERE(node.X,0)
                return node.X,
        Update(self.root)

    def Node_Count(self):
        def Count(node):
            return 1+sum(Count(c) for c in node.children)
        return Count(self.root)

    def Leaf_Count(self):
        def Count(node):
            return (not node.children)+sum(Count(c) for c in node.children)
        return Count(self.root)

    def Volume_Cost(self):
        def Cost(node,parent):
            size=Sphere_Intersection_Volume(node.sphere,parent)
            return size*len(node.children)+sum(Cost(c,node.sphere) for c in node.children)
        return Cost(self.root,self.root.sphere)

    def Sampled_Volume_Cost(self,tolerance=.01):
        def f(X,node=self.root):
            if X not in node.sphere: 
                return 0.
            return len(node.children)+sum(f(X,c) for c in node.children)
        box=self.root.sphere.Bounding_Box()
        return Sobol_Integrate(f,box,tolerance=.01)

    def Optimize(self):
        reduced=[0]
        def Compact(node,parent): 
            node_size=Sphere_Intersection_Volume(node.sphere,parent)
            for c in node.children:
                Compact(c,node.sphere)
            def Child(child):
                if not child.children:
                    return [child]
                child_size=Sphere_Intersection_Volume(child.sphere,node.sphere)
                current_cost=node_size+child_size*len(child.children)+sum(Sphere_Intersection_Volume(c.sphere,child.sphere)*len(c.children) for c in child.children)
                compacted_cost=node_size*len(child.children)+sum(Sphere_Intersection_Volume(c.sphere,node.sphere)*len(c.children) for c in child.children)
                if compacted_cost<=current_cost:
                    reduced[0]=reduced[0]+current_cost-compacted_cost
                    return sum(map(Child,child.children),[])
                else:
                    return [child]
            node.children=sum(map(Child,node.children),[])
        before=self.Volume_Cost()
        Compact(self.root,self.root.sphere)
        after=self.Volume_Cost()
        assert((before-after-reduced[0])/before<.0001)

    def Statistics(self):
        print 'nodes = %r'%self.Node_Count()
        print 'root = %r'%self.root.sphere
        print 'volume cost = %r'%self.Volume_Cost()

    def Print(self):
        def P(node,indent=0):
            print '%s%r, %g'%('  '*indent,node.sphere,node.sphere.Size())
            for c in node.children:
                P(c,indent+1)
        P(self.root)

def Read_Particles(tri_file):
    surface=TRIANGULATED_SURFACE_f.Create()
    Read_From_File("float",tri_file,surface)
    return surface.particles

class TESTS(unittest.TestCase):
    def setUp(self):
        self.rand=RANDOM_NUMBERS_f()
        self.rand.Set_Seed(26732)

    def test_first(self):
        dir=os.path.join(os.environ['PHYSBAM'],'Public_Data/Rigid_Bodies')
        for r in ['box','subdivided_box','plank','Rings_Test/ring_revolve']:
        #for r in ['plank','Rings_Test/ring_revolve']:
            print r
            particles=Read_Particles(os.path.join(dir,r+'.tri'))
            print 'particles = %r'%len(particles)
            hierarchy=SPHERE_HIERARCHY(particles)
            for _ in range(10):
                hierarchy.Incremental_Initialize()
                first_cost=hierarchy.Volume_Cost()
                hierarchy.Optimize()
                optimized_cost=hierarchy.Volume_Cost()
                print 'incremental: count %d, cost %g, optimized cost %g'%(hierarchy.Node_Count(),first_cost,optimized_cost)
            for g in [1,5,10,20]:
                hierarchy.Initialize_From_KD_Tree(g)
                hierarchy.Slow_Update()
                slow_cost=hierarchy.Volume_Cost()
                hierarchy.Optimize()
                optimized_cost=hierarchy.Volume_Cost()
                print 'particles_in_group %d: count %d, slow cost %g, optimized cost %g'%(g,hierarchy.Node_Count(),slow_cost,optimized_cost)
            print

if __name__ == '__main__':
    unittest.main()
