#!/usr/bin/python

from physbam import *

a=LA_i()
for i in range(10): a.Append(i+1)
print a
a[7:8]=[100,200,300]
print a
a[8:200]=-20
print a
a[0:1]=-30
print a
a[4:6]=[12]
print a
