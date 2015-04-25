//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCALE_MATRIX
//##################################################################### 
#ifndef __SCALE_MATRIX__
#define __SCALE_MATRIX__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Matrices/ZERO_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Tools/Vectors/ZERO_VECTOR.h>
#include <cmath>
namespace PhysBAM{
template<class T,int d>
class SCALE_MATRIX
{
public:
    typedef T SCALAR;
    enum WA {m=d,n=d};
    T x;

    explicit SCALE_MATRIX(T z=T()): x(z) {}

    SCALE_MATRIX operator-() const
    {return SCALE_MATRIX(-x);}

    SCALE_MATRIX operator*(T a) const
    {return SCALE_MATRIX(x*a);}

    SCALE_MATRIX operator/(T a) const
    {return SCALE_MATRIX(x/a);}

    SCALE_MATRIX operator+(T a) const
    {return SCALE_MATRIX(x+a);}

    SCALE_MATRIX operator-(T a) const
    {return SCALE_MATRIX(x-a);}

    SCALE_MATRIX Transposed() const
    {return *this;}
};

template<class T,int d> MATRIX<T,d> operator+= (MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s) {return m+=s.x;}
template<class T,int d> MATRIX<T,d> operator-= (MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s) {return m-=s.x;}
template<class T,int d> MATRIX<T,d> operator+= (MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& s) {return m+=1;}
template<class T,int d> MATRIX<T,d> operator-= (MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& s) {return m-=1;}

template<class T,int d> SCALE_MATRIX<T,d> operator+ (const SCALE_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(a.x+b.x);}
template<class T,int d> MATRIX<T,d> operator+ (const SCALE_MATRIX<T,d>& s,const MATRIX<T,d>& m) {return m+s.x;}
template<class T,int d> MATRIX<T,d> operator+ (const MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s) {return m+s.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SCALE_MATRIX<T,d>& s,const SYMMETRIC_MATRIX<T,d>& m) {return m+s.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SYMMETRIC_MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s) {return m+s.x;}
template<class T,int d> SCALE_MATRIX<T,d> operator+ (const ZERO_MATRIX<T,d>& z,const SCALE_MATRIX<T,d>& s) {return s;}
template<class T,int d> SCALE_MATRIX<T,d> operator+ (const SCALE_MATRIX<T,d>& s,const ZERO_MATRIX<T,d>& z) {return s;}

template<class T,int d> SCALE_MATRIX<T,d> operator- (const SCALE_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(a.x-b.x);}
template<class T,int d> MATRIX<T,d> operator- (const SCALE_MATRIX<T,d>& s,const MATRIX<T,d>& m) {return s.x-m;}
template<class T,int d> MATRIX<T,d> operator- (const MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s) {return m-s.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SCALE_MATRIX<T,d>& s,const SYMMETRIC_MATRIX<T,d>& m) {return s.x-m;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SYMMETRIC_MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s) {return m-s.x;}
template<class T,int d> SCALE_MATRIX<T,d> operator- (const ZERO_MATRIX<T,d>& z,const SCALE_MATRIX<T,d>& s) {return -s;}
template<class T,int d> SCALE_MATRIX<T,d> operator- (const SCALE_MATRIX<T,d>& s,const ZERO_MATRIX<T,d>& z) {return s;}

template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Times_Transpose(const ZERO_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Times_Transpose(const SCALE_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}
template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Times_Transpose(const SCALE_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2*a.x*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SCALE_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part()*a.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SCALE_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*(2*a.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a.Twice_Symmetric_Part()*b.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a*(2*b.x);}

template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Transpose_Times(const ZERO_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Transpose_Times(const SCALE_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}
template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Transpose_Times(const SCALE_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2*a.x*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SCALE_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part()*a.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SCALE_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*(2*a.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a.Twice_Symmetric_Part()*b.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a*(2*b.x);}

template<class T,int d> SCALE_MATRIX<T,d> Times_Self_Transpose(const SCALE_MATRIX<T,d>& a) {return a*a;}

template<class T,int d> SCALE_MATRIX<T,d> Transpose_Times_Self(const SCALE_MATRIX<T,d>& a) {return a*a;}

template<class T,int d> VECTOR<T,d> operator*(const SCALE_MATRIX<T,d>& a,const VECTOR<T,d>& b) {return b*a.x;}

template<class T,int m> ZERO_VECTOR<T,m> operator*(const SCALE_MATRIX<T,m>& a,ZERO_VECTOR<T,m> b) {return b;}
template<class T,int m,int n> ZERO_MATRIX<T,m,n> operator*(const ZERO_MATRIX<T,m,n>& a,const SCALE_MATRIX<T,n>& b) {return a;}
template<class T,int m,int n> ZERO_MATRIX<T,m,n> operator*(const SCALE_MATRIX<T,m>& a,const ZERO_MATRIX<T,m,n>& b) {return b;}
template<class T,int d> SCALE_MATRIX<T,d> operator*(const SCALE_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a*b.x;}
template<class T,int m,int n> MATRIX<T,m,n> operator*(const SCALE_MATRIX<T,m>& a,const MATRIX<T,m,n>& b) {return b*a.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const SCALE_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*a.x;}
template<class T,int m,int n> MATRIX<T,m,n> operator*(const MATRIX<T,m,n>& a,const SCALE_MATRIX<T,n>& b) {return a*b.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const SYMMETRIC_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a*b.x;}

}
#endif
