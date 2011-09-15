#!/usr/bin/python
######################################################################
# Copyright 2007, Geoffrey Irving.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# File sobol.py
######################################################################

import unittest
import physbam

RW=physbam.float

class READ_WRITE_TESTS(unittest.TestCase):
    def test_sobol(self):
        m,n=512,512
        image=physbam.ARRAYS_2D_Vf3(1,m,1,n)
        grid=physbam.GRID_Vf2(m,n,-m/2,1.5*m,-n/2,1.5*n,True)
        box=physbam.BOX_Vf3(0,m,0,n,0,m)
        sobol=physbam.SOBOL(box)
        print box

        white=physbam.Vf3(1,1,1)
        black=physbam.Vf3(0,0,0)
        for i in range(1,m+1):
            for j in range(1,n+1):
                image[i,j]=white

        count=m*n/100*10
        #count=10
        for _ in range(count):
            X=sobol.Get_Vector()
            #print X
            assert(box.Thickened(.1).Lazy_Inside(X))
            for i in [1,2,3]:
                cell=grid.Clamp_To_Cell(X.Remove_Index(i))
                image[cell]=image[cell]-1*physbam.Vf3.Axis_Vector(i)

        physbam.IMAGE_f.Write('sobol.png',image)

if __name__ == '__main__':
    unittest.main()
