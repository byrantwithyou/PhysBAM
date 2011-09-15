#!/usr/bin/python

import sys

n=int(sys.argv[1])
t=0
for t in range(n):
    xy=t%4
    z=t/4
    if xy==0 or xy==2: c=1
    else: c=3
    print 'Shader%d="DonutShader%d"'%(t+1,(c-z)%3+1)
