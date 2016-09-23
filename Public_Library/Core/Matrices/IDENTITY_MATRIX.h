//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IDENTITY_MATRIX
//##################################################################### 
#ifndef __IDENTITY_MATRIX__
#define __IDENTITY_MATRIX__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SCALE_MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Core/Vectors/ZERO_VECTOR.h>
#include <cmath>
namespace PhysBAM{
template<class T,int d>
class IDENTITY_MATRIX
{
public:
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

    FIXED_NUMBER<T,d> Trace() const
    {return FIXED_NUMBER<T,d>();}
};

template<class T,int d>
inline std::ostream& operator<< (std::ostream& output_stream,const IDENTITY_MATRIX<T,d>& A)
{return output_stream<<SCALE_MATRIX<T,d>(1);}

template<class T,int m> SCALE_MATRIX<T,m> operator*(const T& a,const IDENTITY_MATRIX<T,m>& b) {return b*a;}

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

template<class T,int d> VECTOR<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const VECTOR<T,d>& b) {return b;}
template<class T,int d> ZERO_VECTOR<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,ZERO_VECTOR<T,d> b) {return b;}

template<class T,int d> SCALE_MATRIX<T,d> operator*(const SCALE_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int d> SCALE_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return b;}
template<class T,int d> ZERO_MATRIX<T,d> operator*(const ZERO_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int d> ZERO_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return b;}
template<class T,int d> IDENTITY_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return b;}
template<class T,int m,int n> MATRIX<T,m,n> operator*(const IDENTITY_MATRIX<T,m>& a,const MATRIX<T,m,n>& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const IDENTITY_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> operator*(const SYMMETRIC_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}
template<class T,int m,int n> MATRIX<T,m,n> operator*(const MATRIX<T,m,n>& a,const IDENTITY_MATRIX<T,n>& b) {return a;}
template<class T,int d> DIAGONAL_MATRIX<T,d> operator*(const DIAGONAL_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a;}

template<class T,int d> T Double_Contract(const SCALE_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a.Trace();}
template<class T,int d> T Double_Contract(const IDENTITY_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b) {return b.Trace();}
template<class T,int d> FIXED_NUMBER<T,0> Double_Contract(const ZERO_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return FIXED_NUMBER<T,0>();}
template<class T,int d> FIXED_NUMBER<T,0> Double_Contract(const IDENTITY_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b) {return FIXED_NUMBER<T,0>();}
template<class T,int d> FIXED_NUMBER<T,d> Double_Contract(const IDENTITY_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return FIXED_NUMBER<T,d>();}
template<class T,int d> T Double_Contract(const IDENTITY_MATRIX<T,d>& a,const MATRIX<T,d>& b) {return b.Trace();}
template<class T,int d> T Double_Contract(const IDENTITY_MATRIX<T,d>& a,const SYMMETRIC_MATRIX<T,d>& b) {return b.Trace();}
template<class T,int d> T Double_Contract(const SYMMETRIC_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a.Trace();}
template<class T,int d> T Double_Contract(const MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b) {return a.Trace();}

template<class T,int d,class OP> auto
Transpose_Times(const IDENTITY_MATRIX<T,d>& a,const OP& b) -> decltype(a*b)
{return a*b;}

template<class T,int m,int n> MATRIX<T,n,m>
Transpose_Times(const MATRIX<T,m,n>& a,const IDENTITY_MATRIX<T,m>& b)
{return a.Transposed();}

template<class T,int m> SYMMETRIC_MATRIX<T,m>
Transpose_Times(const SYMMETRIC_MATRIX<T,m>& a,const IDENTITY_MATRIX<T,m>& b)
{return a;}

template<class T,int m> DIAGONAL_MATRIX<T,m>
Transpose_Times(const DIAGONAL_MATRIX<T,m>& a,const IDENTITY_MATRIX<T,m>& b)
{return a;}

template<class T,int d,class OP> auto
Times_Transpose(const IDENTITY_MATRIX<T,d>& a,const OP& b) -> decltype(b.Transposed())
{return b.Transposed();}

template<class T,int m,int n> MATRIX<T,m,n>
Times_Transpose(const MATRIX<T,m,n>& a,const IDENTITY_MATRIX<T,m>& b)
{return a;}

template<class T,int m> SYMMETRIC_MATRIX<T,m>
Times_Transpose(const SYMMETRIC_MATRIX<T,m>& a,const IDENTITY_MATRIX<T,m>& b)
{return a;}

template<class T,int m> DIAGONAL_MATRIX<T,m>
Times_Transpose(const DIAGONAL_MATRIX<T,m>& a,const IDENTITY_MATRIX<T,m>& b)
{return a;}
}
#endif
