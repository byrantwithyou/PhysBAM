#!/usr/bin/python
import os

examples=[("no_binding",1),("hard_bound",2),("soft_bound_impulse",3),("soft_bound_springs",4)]

for name,num in examples:
    output="%d_%s"%(num,name)
    print "##################################################################################"
    print "Running Example %d for Test %s"%(num,name)
    print "##################################################################################"
    os.system("./solids_3d_nocona -example ARMADILLO_TEST -o cfl-4/%s %d | grep \"END F\""%(output,num))

