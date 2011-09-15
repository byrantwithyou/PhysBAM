#!/usr/bin/python
######################################################################
# Copyright 2007, Geoffrey Irving, Craig Schroeder.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# Class MATRIX_TESTS
######################################################################

import unittest
from physbam import *
from math import *

class MATRIX_TESTS(unittest.TestCase):
    T=float
    epsilon=numeric_limits_f.epsilon()

    def setUp(self):
        self.rand_f=RANDOM_NUMBERS_f()
        self.rand_f.Set_Seed(173281)
        self.rand_d=RANDOM_NUMBERS_d()
        self.rand_d.Set_Seed(173281)
        self.rand=self.rand_f

    ######################################################################
    # Utilities
    ######################################################################

    def Rand(self,T):
        return self.__dict__['rand_'+T.name[0]]

    def Assert_Zero(self,A,tolerance=0,msg=None):
        if isinstance(tolerance,str): tolerance,msg=0,tolerance
        error=A.Max_Abs()
        if not msg and tolerance: msg='error %g, tolerance %g, ratio %g'%(error,tolerance,error/tolerance)
        self.assert_(error<=tolerance,msg)

    def Assert_Equal(self,A,B,tolerance=0,msg=None):
        if isinstance(tolerance,str): tolerance,msg=0,tolerance
        try:
            self.Assert_Zero(A-B,tolerance,msg)
        except AssertionError:
            print 'Assert_Equal failed: %s:' %msg
            print 'A: %r\n%s'%(A,A)
            print 'B: %r\n%s'%(B,B)
            raise

    def Assert_Rotation(self,A,tolerance=None):
        if not tolerance: tolerance=2*A.Rows()*self.epsilon
        try:
            self.Assert_Zero(A.Normal_Equations_Matrix()-1,tolerance)
            if A.Rows()==A.Columns():
                self.assert_(abs(A.Determinant()-1)<=tolerance,'determinant %g'%A.Determinant())
        except AssertionError:
            print 'Assert_Rotation failed:'
            print 'A: %r\n%s'%(A,A)
            print 'A^T A: %r\n%s'%((A.Normal_Equations_Matrix(),)*2)
            raise

    ######################################################################
    # Dynamic tests
    ######################################################################

    def test_dynamic(self):
        size=6
        count=5
        for c in range(count):
            for i in range(1,size+1):
                for j in range(1,size+1):
                    self.Dynamic_Tests_One_Size(i,j)
                    for k in range(1,size+1):
                        self.Dynamic_Tests_Two_Sizes(i,j,k)
                        for m in range(1,size+1):
                            self.Dynamic_Tests_Three_Sizes(i,j,k,m)

    def Dynamic_Tests_One_Size(self,m,n):
        A=MATRIX_MXN_f(m,n)
        B=MATRIX_MXN_f(A)
        D=MATRIX_MXN_f(A)
        self.rand.Fill_Uniform_Matrix(A,-1,1)
        self.rand.Fill_Uniform_Matrix(B,-1,1)
        self.rand.Fill_Uniform_Matrix(D,-1,1)
        s=self.rand.Get_Uniform_Number(-1,1)
        t=self.rand.Get_Uniform_Number(-1,1)
        tolerance=self.epsilon

        self.assert_(A.Rows()==m and A.Columns()==n,"Dimension tests.")
        self.assert_(B.Rows()==m and B.Columns()==n,"Dimension tests.")

        self.Assert_Zero(A-A,"Subtraction with self is zero.")
        self.Assert_Equal(- -A,A,"Negation is its own inverse.")
        self.Assert_Equal(A+A+A,A*3,tolerance*2,"Integer scaling as addition.")
        self.Assert_Equal(-A,A*-1,"Integer scaling as negation.")
        self.Assert_Equal(A*s,s*A,tolerance,"Scalar multiplication commutes.")
        if s: self.Assert_Equal(A*(1/s),A*(1/s),tolerance/abs(s),"Scalar division.")
        self.Assert_Equal(s*(t*A),(s*t)*A,tolerance*2,"Scalar multiplication associates.")

        self.Assert_Equal(A+B,B+A,"Addition commutes.")
        self.Assert_Equal(A-B,A+-B,"Definition of subtraction.")
        self.Assert_Equal(-(A-B),B-A,"Subtraction anticommutes.")
        self.Assert_Equal((B+A)+D,B+(A+D),tolerance*2,"Addition associates.")
        self.Assert_Equal(s*A+s*B,s*(A+B),tolerance*2,"Distributivity of scalar multiplication and addition.")

        self.Assert_Equal(A.Transposed().Transposed(),A,"Transpose is its own inverse.")
        self.Assert_Equal(A.Transposed()+B.Transposed(),(A+B).Transposed(),"Transpose of sum.")
        C=MATRIX_MXN_f(A);C.Transpose();self.Assert_Equal(C,A.Transposed(),"Transpose vs Transposed.")

        C=MATRIX_MXN_f(A);self.Assert_Equal(A,C,"Assignment.")
        C=MATRIX_MXN_f(A);C+=B;self.Assert_Equal(A+B,C,tolerance,"Plus equals.")
        C=MATRIX_MXN_f(A);C-=B;self.Assert_Equal(A-B,C,tolerance,"Minus equals.")
        C=MATRIX_MXN_f(A);C*=s;self.Assert_Equal(A*s,C,tolerance,"Scalar times equals.")
        if s: C=MATRIX_MXN_f(A);C/=s;self.Assert_Equal(A/s,C,tolerance/abs(s),"Scalar divide equals.")

    def Dynamic_Tests_Two_Sizes(self,m,n,p):
        A=MATRIX_MXN_f(m,n)
        B=MATRIX_MXN_f(m,n)
        C=MATRIX_MXN_f(n,p)
        D=MATRIX_MXN_f(n,p)
        self.rand.Fill_Uniform_Matrix(A,-1,1)
        self.rand.Fill_Uniform_Matrix(B,-1,1)
        self.rand.Fill_Uniform_Matrix(C,-1,1)
        self.rand.Fill_Uniform_Matrix(D,-1,1)
        s=self.rand.Get_Uniform_Number(-1,1)
        tolerance=self.epsilon

        self.Assert_Equal((A+B)*C,A*C+B*C,tolerance*A.Columns()*2,"Right distributivity.")
        self.Assert_Equal(A*(C+D),A*C+A*D,tolerance*A.Columns()*2,"Left distributivity (%r %r %r)"%(A,C,D))
        self.Assert_Equal((s*A)*C,s*(A*C),tolerance*A.Columns()*2,"Scalar and matrix multiplication associate.")

        self.Assert_Equal(C.Transposed()*A.Transposed(),(A*C).Transposed(),"Transpose of product.")

    def Dynamic_Tests_Three_Sizes(self,m,n,p,q):
        A=MATRIX_MXN_f(m,n)
        B=MATRIX_MXN_f(n,p)
        C=MATRIX_MXN_f(p,q)
        self.rand.Fill_Uniform_Matrix(A,-1,1)
        self.rand.Fill_Uniform_Matrix(B,-1,1)
        self.rand.Fill_Uniform_Matrix(C,-1,1)
        tolerance=self.epsilon

        self.Assert_Equal((A*B)*C,A*(B*C),tolerance*B.Rows()*B.Columns()*2,"Multiplication associates.")

    ######################################################################
    # Arbitrary (static and dynamic) tests
    ######################################################################

    def test_arbitrary(self):
        def types(m,n):
            t=[MATRIX_MXN_f]
            if m==n or sorted((m,n))==[2,3]:
                t.append(MATRIX(m,n))
            if m==n and m in [2,3]:
                t+=[SYMMETRIC_MATRIX(m),DIAGONAL_MATRIX(m),UPPER_TRIANGULAR_MATRIX(m)]
            return t
        size=3
        count=5
        for i in range(1,size+1):
            for j in range(1,size+1):
                for M1 in types(i,j):
                    for c in range(count):
                        self.Arbitrary_Test_One_Size(M1(INITIAL_SIZE(i),INITIAL_SIZE(j)),count)
                    for k in range(1,size+1):
                        for M2 in types(j,k):
                            for c in range(count):
                                self.Arbitrary_Test_Two_Sizes(M1(INITIAL_SIZE(i),INITIAL_SIZE(j)),M2(INITIAL_SIZE(j),INITIAL_SIZE(k)),count)

    def Conversion_Test(self,A,B):
        if A.__class__.__name__.startswith('MATRIX'):
            self.Assert_Equal(B.__class__(A.__class__(B)),B,"Conversion both ways is exact.")

    def Arbitrary_Test_One_Size(self,A,count):
        M=A.__class__
        B=M(INITIAL_SIZE(A.Rows()),INITIAL_SIZE(A.Columns()))
        C=M(A)
        tolerance=self.epsilon
        self.assert_(B.Columns()==A.Columns() and B.Rows()==A.Rows(),"Dimension tests (B).")
        self.assert_(C.Columns()==A.Columns() and C.Rows()==A.Rows(),"Dimension tests (C).")

        for c in range(count):
            self.rand.Fill_Uniform_Matrix(A,-1,1)
            self.rand.Fill_Uniform_Matrix(B,-1,1)
            D=MATRIX_MXN_f(A)
            E=MATRIX_MXN_f(B)
            s=self.rand.Get_Uniform_Number(-1,1)

            self.Conversion_Test(A,D)

            self.Assert_Equal(MATRIX_MXN_f(-A),-D,tolerance,"Negation matches.")
            self.Assert_Equal(MATRIX_MXN_f(s*A),s*D,tolerance,"Left scaling matches.")
            self.Assert_Equal(MATRIX_MXN_f(A*s),D*s,tolerance,"Right scaling matches.")
            if s: self.Assert_Equal(MATRIX_MXN_f(A/s),D/s,tolerance/abs(s),"Scalar division matches.")
            self.Assert_Equal(MATRIX_MXN_f(A+B),D+E,tolerance,"Addition matches.")
            self.Assert_Equal(MATRIX_MXN_f(A-B),D-E,tolerance,"Subtraction matches.")
            if not M.__name__.startswith('UPPER'):
                self.Assert_Equal(MATRIX_MXN_f(A.Transposed()),D.Transposed(),tolerance,"Transpsoed matches.")

            C=M(A);self.Assert_Equal(MATRIX_MXN_f(C),D,tolerance,"Assignment matches.")
            C=M(A);C*=s;self.Assert_Equal(MATRIX_MXN_f(C),D*s,tolerance,"Scalar times equal matches.")
            if s: C=M(A);C/=s;self.Assert_Equal(MATRIX_MXN_f(C),D/s,tolerance/abs(s),"Scalar divide equal matches.")
            C=M(A);C+=B;self.Assert_Equal(MATRIX_MXN_f(C),D+E,tolerance,"Plus equal matches.")
            C=M(A);C-=B;self.Assert_Equal(MATRIX_MXN_f(C),D-E,tolerance,"Minus equal matches.")

            self.Assert_Equal(MATRIX_MXN_f(A),D,"Inputs not changed.")
            self.Assert_Equal(MATRIX_MXN_f(B),E,"Inputs not changed.")

    def Arbitrary_Test_Two_Sizes(self,A,B,count):
        tolerance=self.epsilon*2
        self.assert_(A.Columns()==B.Rows(),"Dimension tests (A).")

        for c in range(count):
            self.rand.Fill_Uniform_Matrix(A,-1,1)
            self.rand.Fill_Uniform_Matrix(B,-1,1)
            D=MATRIX_MXN_f(A)
            E=MATRIX_MXN_f(B)

            self.Assert_Equal(MATRIX_MXN_f(A*B),D*E,tolerance,"Multiplication matches (%r*%r)."%(A,B))

            self.Assert_Equal(MATRIX_MXN_f(A),D,"Inputs not changed.")
            self.Assert_Equal(MATRIX_MXN_f(B),E,"Inputs not changed.")

    ######################################################################
    # SVD tests
    ######################################################################

    def test_svd(self):
        count=10000
        log_range=-1.5*log(self.epsilon)
        T=float

        def Random_Orthogonal_Matrix(T,m,n):
            Q=self.Rand(T).Get_Rotation(m).Rotation_Matrix()
            if m>n: Q=MATRIX(T,m,n)(Q.Column(1),Q.Column(2))
            return Q

        def Test_Eigensolve(A):
            tolerance=6*self.epsilon
            D,V=A.Fast_Solve_Eigenproblem()
            #print '%r = %r,%r'%(A,D,V)
            self.Assert_Rotation(V)
            self.Assert_Equal(SYMMETRIC_MATRIX(V.Rows()).Conjugate(V,D),A,tolerance*A.Max_Abs())
            #print '%d: %g'%(V.Rows(),(SYMMETRIC_MATRIX(V.Rows()).Conjugate(V,D)-A).Max_Abs()/(tolerance*A.Max_Abs()))

        def Test_SVD(A,U_true=None,D_true=None,V_true=None):
            if min(A.Rows(),A.Columns())==2: tolerance=6*self.epsilon
            else: tolerance=.5*sqrt(self.epsilon)
            U,D,V=A.Fast_Singular_Value_Decomposition()
            try:
                self.Assert_Rotation(U)
                self.Assert_Rotation(V)
                self.Assert_Equal(U*D*V.Transposed(),A,tolerance*A.Max_Abs())
            except AssertionError:
                print '\nSVD Failed:'
                print 'D_true: %r\n%s'%(D_true,D_true)
                print 'D: %r\n%s'%(D,D)
                raise
            #print '%d,%d: %g'%(A.Rows(),A.Columns(),(U*D*V.Transposed()-A).Max_Abs()/(tolerance*A.Max_Abs()))

        for A in [
                MATRIX(T,3,2)(4.656173800800057e-09,-5.5460618277044851e-09,1.3484394045193074e-09,
                              -5.1775682377081315e-09,4.4414298744000848e-09,-4.5764913091182535e-10),
                MATRIX(T,3)(-.46673855799602715,.67466260360310948,.97646986796448998,
                            -.032460753747103721,.046584527749418278,.067431228641151142,
                            -.088885055229687815,.1280389179308779,.18532617511453064),
                MATRIX(T,3)(1e-8,0,0,0,1e-8,0,0,0,1e-8),
                MATRIX(T,3)(2.3588321752040035e-09,-9.6558131480729038e-09,1.0959850449366498e-09,
                            8.8671829608044748e-09,1.6771794267033661e-09,-4.3081475729438225e-09,
                            3.9760504409327016e-09,1.9880497026345722e-09,8.9576046614601957e-09),
                MATRIX(T,3)(.0023588321752040036,-.0096558131480729038,.0010959850449366493,
                            .0088671829608044754,.0016771794267033666,-.0043081475729438235,
                            .003976050440932701,.0019880497026345716,.0089576046614601966),
                MATRIX(T,3)(1,0,0,0,.00036,0,1e-18,0,.00018),
                MATRIX(T,3)(1.3,0,0,0,.0003,0,1e-17,0,0),
                MATRIX(T,3)(1,0,0,0,1e-2,0,0,0,1e-2),
                MATRIX(T,3)(1,0,0,0,1,0,0,0,0),
                MATRIX(T,3)(1,0,0,0,1e-3,0,0,0,1e-6),
                MATRIX(T,3)(1,0,0,0,0,0,0,0,0)]:
            Test_Eigensolve(A.Normal_Equations_Matrix())
            Test_SVD(A)

        progress=PROGRESS_INDICATOR(count)
        for c in range(count):
            progress.Progress()
            for m in [2,3]:
                for n in [2,3]:
                    k=min(m,n)
                    D_true=DIAGONAL_MATRIX(T,k)()
                    for i in range(1,k+1):
                        j=self.Rand(T).Get_Uniform_Integer(0,i-1)
                        if j: D_true[i,i]=D_true[j,j]
                        else: D_true[i,i]=exp(self.Rand(T).Get_Uniform_Number(-log_range,log_range))
                    if self.Rand(T).Get_Number()<.5: D_true[k,k]=-D_true[k,k]
                    V_true=Random_Orthogonal_Matrix(T,n,k)
                    if m==n:
                        A=SYMMETRIC_MATRIX(T,m).Conjugate(V_true,D_true)
                        Test_Eigensolve(A)
                    U_true=Random_Orthogonal_Matrix(T,m,k)
                    A=U_true*D_true*V_true.Transposed()
                    Test_SVD(A,U_true,D_true,V_true)

    ######################################################################
    # Miscellaneous tests
    ######################################################################

    def test_repr(self):
        count=1000
        A=MATRIX(3)()
        for c in range(count):
            self.rand.Fill_Uniform_Matrix(A,-100,100)
            self.assert_(eval(repr(A))==A)

    def test_simplex_minimum_altitude(self):
        A=MATRIX(3)()
        TV=Vf3
        tolerance=10*self.epsilon
        for _ in range(1000):
            self.rand.Fill_Uniform_Matrix(A,-1,1)
            tetrahedron=TETRAHEDRON_f(TV(),A.Column(1),A.Column(2),A.Column(3))
            try:
                assert(abs(A.Simplex_Minimum_Altitude()-tetrahedron.Minimum_Altitude())<tolerance)
            except AssertionError:
                print repr(A)
                print A.Simplex_Minimum_Altitude()
                print tetrahedron.Minimum_Altitude()
                print tolerance
                raise

######################################################################
# Main
######################################################################

if __name__ == '__main__':
    unittest.main()
