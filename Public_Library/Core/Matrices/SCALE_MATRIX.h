//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCALE_MATRIX
//##################################################################### 
#ifndef __SCALE_MATRIX__
#define __SCALE_MATRIX__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Matrices/ZERO_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Core/Vectors/ZERO_VECTOR.h>
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

    SCALE_MATRIX Twice_Symmetric_Part() const
    {return *this*2;}

    SCALE_MATRIX Normal_Equations_Matrix() const
    {return SCALE_MATRIX(x*x);}

    SCALE_MATRIX Outer_Product_Matrix() const
    {return SCALE_MATRIX(x*x);}

    T Trace() const
    {return x*m;}
};

template<class T,int d>
inline std::ostream& operator<< (std::ostream& output_stream,const SCALE_MATRIX<T,d>& A)
{return output_stream<<(SYMMETRIC_MATRIX<T,d>()+A.x);}

template<class T,int m> SCALE_MATRIX<T,m> operator*(const T& a,const SCALE_MATRIX<T,m>& b) {return b*a;}

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
template<class T,int d> DIAGONAL_MATRIX<T,d> operator*(const DIAGONAL_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a*b.x;}

template<class T,int m,int n> FIXED_NUMBER<T,0> Double_Contract(const ZERO_MATRIX<T,m,n>& a,const SCALE_MATRIX<T,n>& b) {return FIXED_NUMBER<T,0>();}
template<class T,int m,int n> FIXED_NUMBER<T,0> Double_Contract(const SCALE_MATRIX<T,m>& a,const ZERO_MATRIX<T,m,n>& b) {return FIXED_NUMBER<T,0>();}
template<class T,int d> T Double_Contract(const SCALE_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return a.x*b.Trace();}
template<class T,int m,int n> T Double_Contract(const SCALE_MATRIX<T,m>& a,const MATRIX<T,m,n>& b) {return a.x*b.Trace();}
template<class T,int d> T Double_Contract(const SCALE_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return a.x*b.Trace();}
template<class T,int m,int n> T Double_Contract(const MATRIX<T,m,n>& a,const SCALE_MATRIX<T,n>& b) {return b.x*a.Trace();}
template<class T,int d> T Double_Contract(const SYMMETRIC_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return b.x*a.Trace();}

template<class T,int d,class OP> auto
Transpose_Times(const SCALE_MATRIX<T,d>& a,const OP& b) -> decltype(a*b)
{return a*b;}

template<class T,int m,int n> MATRIX<T,n,m>
Transpose_Times(const MATRIX<T,m,n>& a,const SCALE_MATRIX<T,m>& b)
{return a.Transposed()*b.x;}

template<class T,int m> SYMMETRIC_MATRIX<T,m>
Transpose_Times(const SYMMETRIC_MATRIX<T,m>& a,const SCALE_MATRIX<T,m>& b)
{return a*b.x;}

template<class T,int m> DIAGONAL_MATRIX<T,m>
Transpose_Times(const DIAGONAL_MATRIX<T,m>& a,const SCALE_MATRIX<T,m>& b)
{return a*b.x;}

template<class T,int d,class OP> auto
Times_Transpose(const SCALE_MATRIX<T,d>& a,const OP& b) -> decltype(a.x*b.Transposed())
{return a.x*b.Transposed();}

template<class T,int m,int n> MATRIX<T,n,m>
Times_Transpose(const MATRIX<T,m,n>& a,const SCALE_MATRIX<T,m>& b)
{return a*b.x;}

template<class T,int m> SYMMETRIC_MATRIX<T,m>
Times_Transpose(const SYMMETRIC_MATRIX<T,m>& a,const SCALE_MATRIX<T,m>& b)
{return a*b.x;}

template<class T,int m> DIAGONAL_MATRIX<T,m>
Times_Transpose(const DIAGONAL_MATRIX<T,m>& a,const SCALE_MATRIX<T,m>& b)
{return a*b.x;}

}
#endif
