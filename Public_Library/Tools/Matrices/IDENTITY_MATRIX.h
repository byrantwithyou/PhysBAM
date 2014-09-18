//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IDENTITY_MATRIX
//##################################################################### 
#ifndef __IDENTITY_MATRIX__
#define __IDENTITY_MATRIX__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SCALE_MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Tools/Vectors/ZERO_VECTOR.h>
#include <cmath>
namespace PhysBAM{
template<class T,int d>
struct IDENTITY_MATRIX
{
    typedef T SCALAR;
    enum WA {m=d,n=d};
    explicit IDENTITY_MATRIX() {}

    SCALE_MATRIX<T,d> operator-() const
    {return SCALE_MATRIX<T,d>(-1);}

    SCALE_MATRIX<T,d> operator*(T a) const
    {return SCALE_MATRIX<T,d>(a);}

    SCALE_MATRIX<T,d> operator/(T a) const
    {return SCALE_MATRIX<T,d>(1/a);}

    IDENTITY_MATRIX Transposed() const
    {return *this;}
};

template<class T,int d> SCALE_MATRIX<T,d> operator+ (const IDENTITY_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2);}
template<class T,int d> SCALE_MATRIX<T,d> operator+ (const IDENTITY_MATRIX<T,d>& i,const SCALE_MATRIX<T,d>& s) {return s+1;}
template<class T,int d> SCALE_MATRIX<T,d> operator+ (const SCALE_MATRIX<T,d>& s,const IDENTITY_MATRIX<T,d>& i) {return s+1;}
template<class T,int d> MATRIX<T,d> operator+ (const IDENTITY_MATRIX<T,d>& i,const MATRIX<T,d>& m) {return m+1;}
template<class T,int d> MATRIX<T,d> operator+ (const MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& i) {return m+1;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const IDENTITY_MATRIX<T,d>& i,const SYMMETRIC_MATRIX<T,d>& m) {return m+1;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SYMMETRIC_MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& i) {return m+1;}
template<class T,int d> IDENTITY_MATRIX<T,d> operator+ (const ZERO_MATRIX<T,d>& z,const IDENTITY_MATRIX<T,d>& s) {return s;}
template<class T,int d> IDENTITY_MATRIX<T,d> operator+ (const IDENTITY_MATRIX<T,d>& s,const ZERO_MATRIX<T,d>& z) {return s;}

template<class T,int d> SCALE_MATRIX<T,d> operator- (const IDENTITY_MATRIX<T,d>& i,const SCALE_MATRIX<T,d>& s) {return SCALE_MATRIX<T,d>(1-s.x);}
template<class T,int d> SCALE_MATRIX<T,d> operator- (const SCALE_MATRIX<T,d>& s,const IDENTITY_MATRIX<T,d>& i) {return s-1;}
template<class T,int d> SCALE_MATRIX<T,d> operator- (const ZERO_MATRIX<T,d>& z,const IDENTITY_MATRIX<T,d>& s) {return SCALE_MATRIX<T,d>(-1);}
template<class T,int d> ZERO_MATRIX<T,d> operator- (const IDENTITY_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return ZERO_MATRIX<T,d>();}
template<class T,int d> MATRIX<T,d> operator- (const IDENTITY_MATRIX<T,d>& i,const MATRIX<T,d>& m) {return (T)1-m;}
template<class T,int d> MATRIX<T,d> operator- (const MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& i) {return m-1;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const IDENTITY_MATRIX<T,d>& i,const SYMMETRIC_MATRIX<T,d>& m) {return (T)1-m;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SYMMETRIC_MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& i) {return m-1;}
template<class T,int d> IDENTITY_MATRIX<T,d> operator- (const IDENTITY_MATRIX<T,d>& s,const ZERO_MATRIX<T,d>& z) {return s;}

template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Times_Transpose(const SCALE_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2*a.x);}
template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Times_Transpose(const IDENTITY_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2);}
template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Times_Transpose(const IDENTITY_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2*b.x);}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Times_Transpose(const IDENTITY_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const IDENTITY_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const IDENTITY_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*2;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Times_Transpose(const ZERO_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a*2;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a.Twice_Symmetric_Part();}

template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Transpose_Times(const ZERO_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Transpose_Times(const IDENTITY_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const IDENTITY_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const IDENTITY_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*2;}
template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Transpose_Times(const SCALE_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2*a.x);}
template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Transpose_Times(const IDENTITY_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2);}
template<class T,int d> SCALE_MATRIX<T,d> Symmetric_Transpose_Times(const IDENTITY_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return SCALE_MATRIX<T,d>(2*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a*2;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a.Twice_Symmetric_Part();}

template<class T,int d> IDENTITY_MATRIX<T,d> Times_Self_Transpose(const IDENTITY_MATRIX<T,d>& a) {return a;}

template<class T,int d> IDENTITY_MATRIX<T,d> Transpose_Times_Self(const IDENTITY_MATRIX<T,d>& a) {return a;}

template<class T,int d> SCALE_MATRIX<T,d> operator*(const SCALE_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int d> SCALE_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return b;}
template<class T,int d> ZERO_MATRIX<T,d> operator*(const ZERO_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int d> VECTOR<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const VECTOR<T,d>& b) {return b;}
template<class T,int d> ZERO_VECTOR<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,ZERO_VECTOR<T,d> b) {return b;}
template<class T,int d> ZERO_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}
template<class T,int d> IDENTITY_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return b;}
template<class T,int d> MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const SYMMETRIC_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int d> MATRIX<T,d> operator*(const MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
}
#endif
