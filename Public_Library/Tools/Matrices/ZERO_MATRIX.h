//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ZERO_MATRIX
//##################################################################### 
#ifndef __ZERO_MATRIX__
#define __ZERO_MATRIX__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Tools/Vectors/ZERO_VECTOR.h>
#include <cmath>
namespace PhysBAM{
template<class T,int mm,int nn>
class ZERO_MATRIX
{
public:
    typedef T SCALAR;
    enum WA {m=mm,n=nn};
    ZERO_MATRIX operator-() const
    {return *this;}

    ZERO_MATRIX operator+(const ZERO_MATRIX&) const
    {return *this;}

    ZERO_MATRIX operator-(const ZERO_MATRIX&) const
    {return *this;}

    ZERO_MATRIX operator*(const T&) const
    {return *this;}

    ZERO_MATRIX operator/(const T&) const
    {return *this;}

    ZERO_MATRIX Transposed() const
    {return *this;}

    template<class OBJ>
    auto Transpose_Times(const OBJ& o) const -> decltype(this->Transposed()*o)
    {return Transposed()*o;}
};

template<class T,int m,int n> ZERO_MATRIX<T,m,n> operator+ (const ZERO_MATRIX<T,m,n>& a,const ZERO_MATRIX<T,m,n>& b) {return a;}
template<class T,int m,int n> MATRIX<T,m,n> operator+ (const ZERO_MATRIX<T,m,n>& z,const MATRIX<T,m,n>& M) {return M;}
template<class T,int m,int n> MATRIX<T,m,n> operator+ (const MATRIX<T,m,n>& M,const ZERO_MATRIX<T,m,n>& z) {return M;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const ZERO_MATRIX<T,d>& z,const SYMMETRIC_MATRIX<T,d>& M) {return M;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator+ (const SYMMETRIC_MATRIX<T,d>& M,const ZERO_MATRIX<T,d>& z) {return M;}

template<class T,int m,int n> ZERO_MATRIX<T,m,n> operator- (const ZERO_MATRIX<T,m,n>& a,const ZERO_MATRIX<T,m,n>& b) {return a;}
template<class T,int m,int n> MATRIX<T,m,n> operator- (const ZERO_MATRIX<T,m,n>& z,const MATRIX<T,m,n>& M) {return -M;}
template<class T,int m,int n> MATRIX<T,m,n> operator- (const MATRIX<T,m,n>& M,const ZERO_MATRIX<T,m,n>& z) {return M;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const ZERO_MATRIX<T,d>& z,const SYMMETRIC_MATRIX<T,d>& M) {return -M;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator- (const SYMMETRIC_MATRIX<T,d>& M,const ZERO_MATRIX<T,d>& z) {return M;}

template<class T,int m,int n> ZERO_MATRIX<T,m,m> Symmetric_Times_Transpose(const ZERO_MATRIX<T,m,n>& a,const ZERO_MATRIX<T,m,n>& b) {return a;}
template<class T,int m,int n> ZERO_MATRIX<T,m,m> Symmetric_Times_Transpose(const ZERO_MATRIX<T,m,n>& a,const MATRIX<T,m,n>& b) {return a;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Times_Transpose(const ZERO_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return a;}
template<class T,int m,int n> ZERO_MATRIX<T,m,m> Symmetric_Times_Transpose(const MATRIX<T,m,n>& a,const ZERO_MATRIX<T,m,n>& b) {return b;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Times_Transpose(const SYMMETRIC_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}

template<class T,int m,int n> ZERO_MATRIX<T,n,n> Symmetric_Transpose_Times(const ZERO_MATRIX<T,m,n>& a,const ZERO_MATRIX<T,m,n>& b) {return a;}
template<class T,int m,int n> ZERO_MATRIX<T,n,n> Symmetric_Transpose_Times(const ZERO_MATRIX<T,m,n>& a,const MATRIX<T,m,n>& b) {return a;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Transpose_Times(const ZERO_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return a;}
template<class T,int m,int n> ZERO_MATRIX<T,n,n> Symmetric_Transpose_Times(const MATRIX<T,m,n>& a,const ZERO_MATRIX<T,m,n>& b) {return b;}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Transpose_Times(const SYMMETRIC_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}

template<class T,int m,int n> ZERO_MATRIX<T,m,m> Times_Self_Transpose(const ZERO_MATRIX<T,m,n>& a) {return a;}

template<class T,int m,int n> ZERO_MATRIX<T,n,n> Transpose_Times_Self(const ZERO_MATRIX<T,m,n>& a) {return a;}
template<class T,int m,int n> SYMMETRIC_MATRIX<T,n> Transpose_Times_Self(const MATRIX<T,m,n>& a) {return a.Normal_Equations_Matrix();}

template<class T,int m,int n> ZERO_VECTOR<T,m> operator*(const ZERO_MATRIX<T,m,n>& a,ZERO_VECTOR<T,n> b) {return b;}
template<class T,int m,int n> ZERO_VECTOR<T,m> operator*(const MATRIX<T,m,n>& a,ZERO_VECTOR<T,n> b) {return b;}
template<class T,int d> ZERO_VECTOR<T,d> operator*(const SYMMETRIC_MATRIX<T,d>& a,ZERO_VECTOR<T,d> b) {return b;}
template<class T,int m,int n> ZERO_VECTOR<T,m> operator*(const ZERO_MATRIX<T,m,n>& a,const VECTOR<T,n>& b) {return ZERO_VECTOR<T,m>();}

template<class T,int m,int n,int p> ZERO_MATRIX<T,m,p> operator*(const ZERO_MATRIX<T,m,n>& a,const ZERO_MATRIX<T,n,p>& b) {return a;}
template<class T,int m,int n,int p> ZERO_MATRIX<T,m,p> operator*(const ZERO_MATRIX<T,m,n>& a,const MATRIX<T,n,p>& b) {return a;}
template<class T,int m,int n> ZERO_MATRIX<T,m,n> operator*(const ZERO_MATRIX<T,m,n>& a,const SYMMETRIC_MATRIX<T,n>& b) {return a;}
template<class T,int m,int n,int p> ZERO_MATRIX<T,m,p> operator*(const MATRIX<T,m,n>& a,const ZERO_MATRIX<T,n,p>& b) {return b;}
template<class T,int m,int n> ZERO_MATRIX<T,m,n> operator*(const SYMMETRIC_MATRIX<T,m>& a,const ZERO_MATRIX<T,m,n>& b) {return b;}

template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Outer_Product(ZERO_VECTOR<T,d> u,const VECTOR<T,d>& v) {return ZERO_MATRIX<T,d>();}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Outer_Product(ZERO_VECTOR<T,d> u,ZERO_VECTOR<T,d> v) {return ZERO_MATRIX<T,d>();}
template<class T,int d> ZERO_MATRIX<T,d> Symmetric_Outer_Product(const VECTOR<T,d>& u,ZERO_VECTOR<T,d> v) {return ZERO_MATRIX<T,d>();}

template<class T,int d> ZERO_MATRIX<T,d> Outer_Product(const ZERO_VECTOR<T,d>& u) {return ZERO_MATRIX<T,d>();}
template<class T,int d> ZERO_MATRIX<T,d> Outer_Product(const ZERO_VECTOR<T,d>& u,const VECTOR<T,d>& v) {return ZERO_MATRIX<T,d>();}
template<class T,int d> ZERO_MATRIX<T,d> Outer_Product(const VECTOR<T,d>& u,const ZERO_VECTOR<T,d>& v) {return ZERO_MATRIX<T,d>();}
template<class T,int d> ZERO_MATRIX<T,d> Outer_Product(const ZERO_VECTOR<T,d>& u,const ZERO_VECTOR<T,d>& v) {return ZERO_MATRIX<T,d>();}

template<class T,int d,class OP> auto
Transpose_Times(const ZERO_MATRIX<T,d>& a,const OP& b) -> decltype(a*b)
{return a*b;}

template<class T,int m,int n,int p> ZERO_MATRIX<T,n,p>
Transpose_Times(const MATRIX<T,m,n>& a,const ZERO_MATRIX<T,m,p>& b)
{return ZERO_MATRIX<T,n,p>();}

template<class T,int m,int n> ZERO_MATRIX<T,m,n>
Transpose_Times(const SYMMETRIC_MATRIX<T,m>& a,const ZERO_MATRIX<T,m,n>& b)
{return ZERO_MATRIX<T,m,n>();}

}
#endif
