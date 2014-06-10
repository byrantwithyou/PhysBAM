//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRIMITIVE_MATRICES
//##################################################################### 
#ifndef __PRIMITIVE_MATRICES__
#define __PRIMITIVE_MATRICES__

#include <Tools/Auto_Diff/PRIMITIVE_VECTORS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV>
struct ZERO_MAT
{
    typedef typename TV::SCALAR SCALAR;typedef SCALAR T;
    enum WA {m=TV::m,n=TV::m};
    ZERO_MAT operator-() const
    {return *this;}

    ZERO_MAT operator*(const T&) const
    {return *this;}

    ZERO_MAT operator/(const T&) const
    {return *this;}

    ZERO_MAT Transposed() const
    {return *this;}
};

template<class TV>
struct SCALE_MAT
{
    typedef typename TV::SCALAR SCALAR;typedef SCALAR T;
    enum WA {m=TV::m,n=TV::m};
    T x;

    explicit SCALE_MAT(T z=T()): x(z) {}

    SCALE_MAT operator-() const
    {return SCALE_MAT(-x);}

    SCALE_MAT operator*(T a) const
    {return SCALE_MAT(x*a);}

    SCALE_MAT operator/(T a) const
    {return SCALE_MAT(x/a);}

    SCALE_MAT operator+(T a) const
    {return SCALE_MAT(x+a);}

    SCALE_MAT operator-(T a) const
    {return SCALE_MAT(x-a);}

    SCALE_MAT Transposed() const
    {return *this;}
};

template<class TV>
struct ID_MAT
{
    typedef typename TV::SCALAR SCALAR;typedef SCALAR T;
    enum WA {m=TV::m,n=TV::m};
    explicit ID_MAT() {}

    SCALE_MAT<TV> operator-() const
    {return SCALE_MAT<TV>(-1);}

    SCALE_MAT<TV> operator*(T a) const
    {return SCALE_MAT<TV>(a);}

    SCALE_MAT<TV> operator/(T a) const
    {return SCALE_MAT<TV>(1/a);}

    ID_MAT Transposed() const
    {return *this;}
};

template<class T> struct IS_MATRIX{static const int value=0;};
template<class TV> struct IS_MATRIX<ZERO_MAT<TV> > {static const int value=1;};
template<class T,int d> struct IS_MATRIX<MATRIX<T,d> > {static const int value=1;};
template<class T,int d> struct IS_MATRIX<SYMMETRIC_MATRIX<T,d> > {static const int value=1;};
template<class TV> struct IS_MATRIX<ID_MAT<TV> > {static const int value=1;};
template<class TV> struct IS_MATRIX<SCALE_MAT<TV> > {static const int value=1;};

template<class T> struct IS_SYM_MATRIX{static const int value=0;};
template<class TV> struct IS_SYM_MATRIX<ZERO_MAT<TV> > {static const int value=1;};
template<class T,int d> struct IS_SYM_MATRIX<SYMMETRIC_MATRIX<T,d> > {static const int value=1;};
template<class TV> struct IS_SYM_MATRIX<ID_MAT<TV> > {static const int value=1;};
template<class TV> struct IS_SYM_MATRIX<SCALE_MAT<TV> > {static const int value=1;};

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator+ (const ZERO_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const SCALE_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(a.x+b.x);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const ID_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2);}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator- (const ZERO_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const SCALE_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(a.x-b.x);}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator- (const ID_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return ZERO_MAT<VECTOR<T,d> >();}

template<class T,int d> MATRIX<T,d> operator+ (const ZERO_MAT<VECTOR<T,d> >& z,const MATRIX<T,d>& m) {return m;}
template<class T,int d> MATRIX<T,d> operator+ (const MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z) {return m;}
template<class T,int d> MATRIX<T,d> operator- (const ZERO_MAT<VECTOR<T,d> >& z,const MATRIX<T,d>& m) {return -m;}
template<class T,int d> MATRIX<T,d> operator- (const MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z) {return m;}

template<class T,int d> MATRIX<T,d> operator+ (const ID_MAT<VECTOR<T,d> >& i,const MATRIX<T,d>& m) {return m+1;}
template<class T,int d> MATRIX<T,d> operator+ (const MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i) {return m+1;}
template<class T,int d> MATRIX<T,d> operator- (const ID_MAT<VECTOR<T,d> >& i,const MATRIX<T,d>& m) {return (T)1-m;}
template<class T,int d> MATRIX<T,d> operator- (const MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i) {return m-1;}

template<class T,int d> MATRIX<T,d> operator+ (const SCALE_MAT<VECTOR<T,d> >& s,const MATRIX<T,d>& m) {return m+s.x;}
template<class T,int d> MATRIX<T,d> operator+ (const MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s) {return m+s.x;}
template<class T,int d> MATRIX<T,d> operator- (const SCALE_MAT<VECTOR<T,d> >& s,const MATRIX<T,d>& m) {return s.x-m;}
template<class T,int d> MATRIX<T,d> operator- (const MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s) {return m-s.x;}

template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const ZERO_MAT<VECTOR<T,d> >& z,const SYMMETRIC_MATRIX<T,d>& m) {return m;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SYMMETRIC_MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z) {return m;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const ZERO_MAT<VECTOR<T,d> >& z,const SYMMETRIC_MATRIX<T,d>& m) {return -m;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SYMMETRIC_MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z) {return m;}

template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const ID_MAT<VECTOR<T,d> >& i,const SYMMETRIC_MATRIX<T,d>& m) {return m+1;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SYMMETRIC_MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i) {return m+1;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const ID_MAT<VECTOR<T,d> >& i,const SYMMETRIC_MATRIX<T,d>& m) {return (T)1-m;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SYMMETRIC_MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i) {return m-1;}

template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SCALE_MAT<VECTOR<T,d> >& s,const SYMMETRIC_MATRIX<T,d>& m) {return m+s.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SYMMETRIC_MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s) {return m+s.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SCALE_MAT<VECTOR<T,d> >& s,const SYMMETRIC_MATRIX<T,d>& m) {return s.x-m;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SYMMETRIC_MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s) {return m-s.x;}

template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const ZERO_MAT<VECTOR<T,d> >& z,const SCALE_MAT<VECTOR<T,d> >& s) {return s;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const SCALE_MAT<VECTOR<T,d> >& s,const ZERO_MAT<VECTOR<T,d> >& z) {return s;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const ZERO_MAT<VECTOR<T,d> >& z,const SCALE_MAT<VECTOR<T,d> >& s) {return -s;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const SCALE_MAT<VECTOR<T,d> >& s,const ZERO_MAT<VECTOR<T,d> >& z) {return s;}

template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const ID_MAT<VECTOR<T,d> >& i,const SCALE_MAT<VECTOR<T,d> >& s) {return s+1;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const SCALE_MAT<VECTOR<T,d> >& s,const ID_MAT<VECTOR<T,d> >& i) {return s+1;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const ID_MAT<VECTOR<T,d> >& i,const SCALE_MAT<VECTOR<T,d> >& s) {return SCALE_MAT<VECTOR<T,d> >(1-s.x);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const SCALE_MAT<VECTOR<T,d> >& s,const ID_MAT<VECTOR<T,d> >& i) {return s-1;}

template<class T,int d> ID_MAT<VECTOR<T,d> > operator+ (const ZERO_MAT<VECTOR<T,d> >& z,const ID_MAT<VECTOR<T,d> >& s) {return s;}
template<class T,int d> ID_MAT<VECTOR<T,d> > operator+ (const ID_MAT<VECTOR<T,d> >& s,const ZERO_MAT<VECTOR<T,d> >& z) {return s;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const ZERO_MAT<VECTOR<T,d> >& z,const ID_MAT<VECTOR<T,d> >& s) {return SCALE_MAT<VECTOR<T,d> >(-1);}
template<class T,int d> ID_MAT<VECTOR<T,d> > operator- (const ID_MAT<VECTOR<T,d> >& s,const ZERO_MAT<VECTOR<T,d> >& z) {return s;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ZERO_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ZERO_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ZERO_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ZERO_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ZERO_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return a;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ID_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ID_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const ID_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const ID_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const ID_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*2;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const SCALE_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const SCALE_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2*a.x);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const SCALE_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2*a.x*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SCALE_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part()*a.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SCALE_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*(2*a.x);}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const MATRIX<T,d>& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const MATRIX<T,d>& a,const ID_MAT<VECTOR<T,d> >& b) {return a.Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const MATRIX<T,d>& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a.Twice_Symmetric_Part()*b.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const MATRIX<T,d>& a,const MATRIX<T,d>& b) {return a.Times_Transpose(b).Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return a.Times_Transpose(b).Twice_Symmetric_Part();}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const ID_MAT<VECTOR<T,d> >& b) {return a*2;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a*(2*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return a.Times_Transpose(b).Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return a.Times_Transpose(b).Twice_Symmetric_Part();}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ZERO_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ZERO_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ZERO_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ZERO_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ZERO_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return a;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ID_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ID_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const ID_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const ID_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const ID_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*2;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const SCALE_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const SCALE_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2*a.x);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const SCALE_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2*a.x*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SCALE_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return b.Twice_Symmetric_Part()*a.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SCALE_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*(2*a.x);}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const MATRIX<T,d>& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const MATRIX<T,d>& a,const ID_MAT<VECTOR<T,d> >& b) {return a.Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const MATRIX<T,d>& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a.Twice_Symmetric_Part()*b.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const MATRIX<T,d>& a,const MATRIX<T,d>& b) {return a.Transpose_Times(b).Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return a.Transpose_Times(b).Twice_Symmetric_Part();}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const ID_MAT<VECTOR<T,d> >& b) {return a*2;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a*(2*b.x);}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return a.Transpose_Times(b).Twice_Symmetric_Part();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return a.Transpose_Times(b).Twice_Symmetric_Part();}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Times_Self_Transpose(const ZERO_MAT<VECTOR<T,d> >& a) {return a;}
template<class T,int d> ID_MAT<VECTOR<T,d> > Times_Self_Transpose(const ID_MAT<VECTOR<T,d> >& a) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Times_Self_Transpose(const SCALE_MAT<VECTOR<T,d> >& a) {return a*a;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Times_Self_Transpose(const MATRIX<T,d>& a) {return a.Outer_Product_Matrix();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Times_Self_Transpose(const SYMMETRIC_MATRIX<T,d>& a) {return a*a;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Transpose_Times_Self(const ZERO_MAT<VECTOR<T,d> >& a) {return a;}
template<class T,int d> ID_MAT<VECTOR<T,d> > Transpose_Times_Self(const ID_MAT<VECTOR<T,d> >& a) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Transpose_Times_Self(const SCALE_MAT<VECTOR<T,d> >& a) {return a*a;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Transpose_Times_Self(const MATRIX<T,d>& a) {return a.Normal_Equations_Matrix();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Transpose_Times_Self(const SYMMETRIC_MATRIX<T,d>& a) {return a*a;}

template<class T,int d> ZERO_VEC<VECTOR<T,d> > operator*(const ZERO_MAT<VECTOR<T,d> >& a,ZERO_VEC<VECTOR<T,d> > b) {return b;}
template<class T,int d> ZERO_VEC<VECTOR<T,d> > operator*(const ID_MAT<VECTOR<T,d> >& a,ZERO_VEC<VECTOR<T,d> > b) {return b;}
template<class T,int d> ZERO_VEC<VECTOR<T,d> > operator*(const SCALE_MAT<VECTOR<T,d> >& a,ZERO_VEC<VECTOR<T,d> > b) {return b;}
template<class T,int d> ZERO_VEC<VECTOR<T,d> > operator*(const MATRIX<T,d>& a,ZERO_VEC<VECTOR<T,d> > b) {return b;}
template<class T,int d> ZERO_VEC<VECTOR<T,d> > operator*(const SYMMETRIC_MATRIX<T,d>& a,ZERO_VEC<VECTOR<T,d> > b) {return b;}

template<class T,int d> ZERO_VEC<VECTOR<T,d> > operator*(const ZERO_MAT<VECTOR<T,d> >& a,const VECTOR<T,d>& b) {return ZERO_VEC<VECTOR<T,d> >();}
template<class T,int d> VECTOR<T,d> operator*(const ID_MAT<VECTOR<T,d> >& a,const VECTOR<T,d>& b) {return b;}
template<class T,int d> VECTOR<T,d> operator*(const SCALE_MAT<VECTOR<T,d> >& a,const VECTOR<T,d>& b) {return b*a.x;}
template<class T,int d> VECTOR<T,d> operator*(const MATRIX<T,d>& a,const VECTOR<T,d>& b) {return a*b;}
template<class T,int d> VECTOR<T,d> operator*(const SYMMETRIC_MATRIX<T,d>& a,const VECTOR<T,d>& b) {return a*b;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const ZERO_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const ZERO_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const ZERO_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const ZERO_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return a;}
template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const ZERO_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return a;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const ID_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> ID_MAT<VECTOR<T,d> > operator*(const ID_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator*(const ID_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> MATRIX<T,d> operator*(const ID_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const ID_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return b;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const SCALE_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator*(const SCALE_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator*(const SCALE_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a*b.x;}
template<class T,int d> MATRIX<T,d> operator*(const SCALE_MAT<VECTOR<T,d> >& a,const MATRIX<T,d>& b) {return b*a.x;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const SCALE_MAT<VECTOR<T,d> >& a,const SYMMETRIC_MATRIX<T,d>& b) {return b*a.x;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const MATRIX<T,d>& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> MATRIX<T,d> operator*(const MATRIX<T,d>& a,const ID_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> MATRIX<T,d> operator*(const MATRIX<T,d>& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a*b.x;}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator*(const SYMMETRIC_MATRIX<T,d>& a,const ZERO_MAT<VECTOR<T,d> >& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const SYMMETRIC_MATRIX<T,d>& a,const ID_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const SYMMETRIC_MATRIX<T,d>& a,const SCALE_MAT<VECTOR<T,d> >& b) {return a*b.x;}

template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(ZERO_MAT<VECTOR<T,d> > z){return SYMMETRIC_MATRIX<T,d>();}
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(SCALE_MAT<VECTOR<T,d> > z){return SYMMETRIC_MATRIX<T,d>()+z.x;}
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(ID_MAT<VECTOR<T,d> > z){return SYMMETRIC_MATRIX<T,d>::Identity_Matrix();}
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(SYMMETRIC_MATRIX<T,d> z){return z;}
template<class T,int d> inline MATRIX<T,d> Cast_Helper(const MATRIX<T,d>& z){return z;}

template<class T> typename ENABLE_IF<IS_MATRIX<T>::value>::TYPE Fill_From(T& a,const T& b){a=b;}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const SYMMETRIC_MATRIX<T,d>& z){m=z;}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z){m=MATRIX<T,d>();}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i){m=MATRIX<T,d>::Identity_Matrix();}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s){m=MATRIX<T,d>()+s.x;}
template<class T,int d> void Fill_From(SYMMETRIC_MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z){m=SYMMETRIC_MATRIX<T,d>();}
template<class T,int d> void Fill_From(SYMMETRIC_MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i){m=SYMMETRIC_MATRIX<T,d>::Identity_Matrix();}
template<class T,int d> void Fill_From(SYMMETRIC_MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s){m=SYMMETRIC_MATRIX<T,d>()+s.x;}
template<class T,int d> void Fill_From(SCALE_MAT<VECTOR<T,d> >& s,const ZERO_MAT<VECTOR<T,d> >& z) {s.x=0;}
template<class T,int d> void Fill_From(SCALE_MAT<VECTOR<T,d> >& s,const ID_MAT<VECTOR<T,d> >& i) {s.x=1;}


template<class T,int d> SYMMETRIC_MATRIX<T,d>
Symmetric_Outer_Product_Helper(const VECTOR<T,d>& u,const VECTOR<T,d>& v)
{return SYMMETRIC_MATRIX<T,d>::Symmetric_Outer_Product(u,v);}

template<class TV> ZERO_MAT<TV> Symmetric_Outer_Product_Helper(ZERO_VEC<TV> u,const TV& v) {return ZERO_MAT<TV>();}
template<class TV> ZERO_MAT<TV> Symmetric_Outer_Product_Helper(ZERO_VEC<TV> u,ZERO_VEC<TV> v) {return ZERO_MAT<TV>();}
template<class TV> ZERO_MAT<TV> Symmetric_Outer_Product_Helper(const TV& u,ZERO_VEC<TV> v) {return ZERO_MAT<TV>();}

template<class T,int d> void
Outer_Product_Helper(SYMMETRIC_MATRIX<T,d>& r,const VECTOR<T,d>& u)
{r=SYMMETRIC_MATRIX<T,d>::Outer_Product(u);}

template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Outer_Product_Helper(const TV& u) {return SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(u);}
template<class TV> ZERO_MAT<TV> Outer_Product_Helper(const ZERO_VEC<TV>& u) {return ZERO_MAT<TV>();}

template<class TV> MATRIX<typename TV::SCALAR,TV::m> Outer_Product_Helper(const TV& u,const TV& v) {return MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(u,v);}
template<class TV> ZERO_MAT<TV> Outer_Product_Helper(const ZERO_VEC<TV>& u,const TV& v) {return ZERO_MAT<TV>();}
template<class TV> ZERO_MAT<TV> Outer_Product_Helper(const TV& u,const ZERO_VEC<TV>& v) {return ZERO_MAT<TV>();}
template<class TV> ZERO_MAT<TV> Outer_Product_Helper(const ZERO_VEC<TV>& u,const ZERO_VEC<TV>& v) {return ZERO_MAT<TV>();}

template<class T_MAT> typename ENABLE_IF<IS_MATRIX<T_MAT>::value,T_MAT>::TYPE Choose(const T_MAT& a,const T_MAT& b);
template<class T_MAT,class T_MAT1> typename ENABLE_IF<IS_MATRIX<T_MAT>::value&&IS_MATRIX<T_MAT1>::value&&(!IS_SYM_MATRIX<T_MAT>::value || !IS_SYM_MATRIX<T_MAT1>::value),MATRIX<typename T_MAT::SCALAR,T_MAT::m> >::TYPE Choose(const T_MAT& a,const T_MAT1& b);
template<class T_MAT,class T_MAT1> typename ENABLE_IF<IS_SYM_MATRIX<T_MAT>::value && IS_SYM_MATRIX<T_MAT1>::value,SYMMETRIC_MATRIX<typename T_MAT::SCALAR,T_MAT::m> >::TYPE Choose(const T_MAT& a,const T_MAT1& b);
template<class TV,class T_MAT> T_MAT Choose(const T_MAT& a,const ZERO_MAT<TV>& b);
template<class TV,class T_MAT> T_MAT Choose(const ZERO_MAT<TV>& a,const T_MAT& b);
template<class TV> ZERO_MAT<TV> Choose(const ZERO_MAT<TV>& a,const ZERO_MAT<TV>& b);
template<class TV> SCALE_MAT<TV> Choose(const SCALE_MAT<TV>& a,const ID_MAT<TV>& b);
template<class TV> SCALE_MAT<TV> Choose(const ID_MAT<TV>& a,const SCALE_MAT<TV>& b);
}
}
#endif
