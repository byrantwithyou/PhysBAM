#!/usr/bin/python
######################################################################
# Copyright 2008, Geoffrey Irving.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# Scripts particles.py
######################################################################

import unittest
import tempfile
import shutil
import copy

from physbam import *
from math import *

class particles_tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.mkdtemp(prefix='particles')

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_read_write(self):
        file=self.tmp+'/particles.gz'
        # create particles
        particles=SOLIDS_PARTICLE_Vf3()
        particles.Add_Particles(10)
        particles.Add_Attribute(PARTICLE_VELOCITY_ATTRIBUTE_Vf3())
        particles.Add_Attribute(PARTICLE_RADIUS_ATTRIBUTE_f())
        for i in range(1,10+1):
            particles.X[i].x=i
            particles.V[i].y=-i
            particles.radius[i]=2**i
        # write particles
        Write_To_File(float,file,particles)
        # read particles
        new_particles=SOLIDS_PARTICLE_Vf3()
        Read_From_File(float,file,new_particles)
        # compare particles
        assert particles.X==new_particles.X
        assert particles.V==new_particles.V
        for i in range(1,10+1):
            assert particles.radius[i]==2**i
        assert particles==new_particles

    def test_print(self):
        particles=PARTICLE_LEVELSET_PARTICLE_Vf2()
        particles.Add_Attribute(PARTICLE_ENERGY_ATTRIBUTE_f())
        particles.Add_Particle()
        particles.X[1]=(1,2)
        particles.E[1]=123
        assert len(particles)==1
        assert str(particles[1])

    def test_add(self):
        particles=SOLIDS_PARTICLE_Vf3()
        def expect(size):
            assert len(particles)==len(particles.X)==size
        for i in range(3): 
            expect(i)
            particles.Add_Particle()
        particles.Add_Particles(2)
        expect(5)
        del particles[2]
        expect(4)

    def test_attributes(self):
        first=SOLIDS_PARTICLE_Vf3()
        second=SOLIDS_PARTICLE_Vf3()
        second.Add_Attribute(PARTICLE_RADIUS_ATTRIBUTE_f())
        assert first.Number_Of_Attributes()==1
        assert second.Number_Of_Attributes()==2
        first.Add_Attributes(second)
        assert first.Number_Of_Attributes()==2
        first.Add_Attributes(second)
        assert first.Number_Of_Attributes()==2

    def test_clone(self):
        particles=SOLIDS_PARTICLE_Vf3()
        particles.Add_Attribute(PARTICLE_RADIUS_ATTRIBUTE_f())
        particles.Add_Particle()
        particles.radius[1]=3
        clone=copy.copy(particles)
        assert particles==clone
        particles.Add_Particle()
        assert particles!=clone

    def test_copy(self):
        TV=Vf3
        normal=PARTICLE_LEVELSET_PARTICLE_Vf3()
        removed=PARTICLE_LEVELSET_REMOVED_PARTICLE_Vf3()
        normal.Add_Particle()
        normal.Add_Attribute(PARTICLE_ENERGY_ATTRIBUTE_f())
        normal.Add_Attribute(PARTICLE_DENSITY_ATTRIBUTE_f())
        removed.Add_Particle()
        removed.Add_Attribute(PARTICLE_DENSITY_ATTRIBUTE_f())
        normal.X[1]=(1,2,3)
        normal.rho[1]=1000
        removed.V[1]=(4,5,6)
        removed[1]=normal[1]
        assert removed.X[1]==(1,2,3)
        assert removed.rho[1]==1000
        assert removed.V[1]==TV()
        removed.X[1]=(3,2,1)
        removed.rho[1]=2000
        normal.E[1]=10
        normal[1]=removed[1]
        assert normal.X[1]==(3,2,1)
        assert normal.rho[1]==2000
        assert not normal.E[1]

######################################################################
# Main
######################################################################

if __name__ == '__main__':
    unittest.main()
