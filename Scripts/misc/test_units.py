#!/usr/bin/python

# ./test_units.py <program to run with arguments>

import os
import subprocess
import re
import sys
import math

tol=1e-10

class vec:
    def __init__(self,a):
        self.x=a
    def __add__(self,v):
        return vec([i[0]+i[1] for i in zip(self.x,v.x)])
    def __sub__(self,v):
        return vec([i[0]-i[1] for i in zip(self.x,v.x)])
    def __mul__(self,a):
        return vec([a*i for i in self.x])
    def dot(self,v):
        return sum([i[0]*i[1] for i in zip(self.x,v.x)])
    def norm(self):
        return self.dot(self)

tu=re.compile("TEST_UNITS (\S+) (\S+) (.*?) @BEGIN@ (.*?) @END@")
nu=re.compile("[0-9eE.+-]+")
def parse_sim(m,k,s):
    print("EXEC %g %g %g"%(m,k,s));
    s=subprocess.check_output(sys.argv[1:]+['-m',str(m),'-kg',str(k),'-s',str(s)])
    r = []
    for g in tu.findall(s):
        r.append((g[0],g[1],g[2],vec([float(i) for i in nu.findall(g[3])])))
    return r

X=parse_sim(1,1,1)
M=parse_sim(2,1,1)
K=parse_sim(1,2,1)
S=parse_sim(1,1,2)
print "EXEC DONE"
N=min(len(X),len(M),len(K),len(S))

def get_pow(x,u,name):
    r=u.dot(x)/x.norm()
    e=(u-x*r).norm()/u.norm()
    fail=None
    if e>tol*tol:
        print "Unit error (%s): %g"%(name,e)
        return None
    p=math.log(r)/math.log(2)
    if abs(p-round(p))>tol:
        print "Power (%s) not a unit: %g"%(name,p)
    return int(round(p))

for i in range(N):
    x=X[i][3]
    m=M[i][3]
    k=K[i][3]
    s=S[i][3]
    if(len(x.x)!=len(m.x) or len(x.x)!=len(k.x) or len(x.x)!=len(s.x)):
        print "Different numbers of entries\nfile %s\nfunction %s\nquantity %s"%X[i][0:3]
        exit()
    if x.norm()<tol*tol:
        print "%s: 0"%X[i][2]
        continue
    pm=get_pow(x,m,'m')
    pk=get_pow(x,k,'kg')
    ps=get_pow(x,s,'s')
    if pm==None or pk==None or ps==None:
        print "No clean unit\nfile %s\nfunction %s\nquantity %s"%X[i][0:3]
        print "Got:",pm,pk,ps
        exit()
    seen=None
    sys.stdout.write(X[i][2]+":")
    if pm==0 and pk==0 and ps==0:
        print " 1"
    else:
        if pk==1: sys.stdout.write(" kg")
        elif pk!=0: sys.stdout.write(" kg^%i"%pk)
        if pm==1: sys.stdout.write(" m")
        elif pm!=0: sys.stdout.write(" m^%i"%pm)
        if ps==1: sys.stdout.write(" s")
        elif ps!=0: sys.stdout.write(" s^%i"%ps)
        print

if len(X)!=len(M) or len(X)!=len(K) or len(X)!=len(S):
    print "different number of lines: %i %i %i %i\n"%(len(X),len(M),len(K),len(S))
