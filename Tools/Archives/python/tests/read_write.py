#!/usr/bin/python
######################################################################
# Copyright 2007, Geoffrey Irving.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# File read_write.py
######################################################################

import unittest
import tempfile
import physbam

RW=physbam.float

class READ_WRITE_TESTS(unittest.TestCase):
    def setUp(self):
        self.file=tempfile.NamedTemporaryFile()

    def test_read_write(self):
        objects=[physbam.LA_LA_i([[1,2,3],[4,5,6]])]
        for a in objects:
            physbam.Write_To_File(RW,self.file.name,a)
            b=type(a)()
            physbam.Read_From_File(RW,self.file.name,b)
            assert(a==b)

if __name__ == '__main__':
    unittest.main()
