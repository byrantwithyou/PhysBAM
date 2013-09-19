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

template<int m,int s> struct M_FLAGS;
template<class TV>
struct ZERO_MAT
{
    typedef M_FLAGS<0,0> FLAGS;
    typedef typename TV::SCALAR T;
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
    typedef M_FLAGS<1,1> FLAGS;
    typedef typename TV::SCALAR T;
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
    typedef M_FLAGS<0,1> FLAGS;
    typedef typename TV::SCALAR T;
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

template<int m,int s> struct M_FLAGS;
template<int nn,int m,int s> struct CM_FLAGS;
template<int nn,int m,int s> struct MM_FLAGS;

template<class A,class B> struct H_ADD
{
    enum WA {mask=A::mask|B::mask|(A::sp_mask&B::sp_mask),sp_mask=(~A::mask|A::sp_mask)&(~B::mask|B::sp_mask)&(A::mask|B::mask|A::sp_mask|B::sp_mask)};
    typedef CM_FLAGS<A::n,mask,sp_mask> CM;
    typedef MM_FLAGS<A::n,mask,sp_mask> MM;
    typedef M_FLAGS<mask,sp_mask> M;
};

template<int m,int s> struct M_FLAGS
{
    enum WA{mask=m,sp_mask=s};

    M_FLAGS<m,s> operator+() const {return M_FLAGS<m,s>();}
    M_FLAGS<m|s,s> operator-() const {return M_FLAGS<m|s,s>();}

    template<int m2,int s2> typename H_ADD<M_FLAGS,M_FLAGS<m2,s2> >::M operator+(M_FLAGS<m2,s2> x) const
    {return typename H_ADD<M_FLAGS,M_FLAGS<m2,s2> >::M();}
    template<int m2,int s2> typename H_ADD<M_FLAGS,M_FLAGS<s2|m2,s2> >::M operator-(M_FLAGS<m2,s2> x) const
    {return typename H_ADD<M_FLAGS,M_FLAGS<s2|m2,s2> >::M();}
    M_FLAGS<m|s,s> operator*(double a) const {return M_FLAGS<m|s,s>();}
    M_FLAGS<m|s,s> operator/(double a) const {return M_FLAGS<m|s,s>();}
};
template<class T,int d> struct GET_FLAGS<MATRIX<T,d> > {typedef M_FLAGS<1,0> TYPE;};
template<class T,int d> struct GET_FLAGS<SYMMETRIC_MATRIX<T,d> > {typedef M_FLAGS<1,0> TYPE;};

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator+ (const ZERO_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const SCALE_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(a.x+b.x);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator+ (const ID_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(2);}

template<class T,int d> ZERO_MAT<VECTOR<T,d> > operator- (const ZERO_MAT<VECTOR<T,d> >& a,const ZERO_MAT<VECTOR<T,d> >& b) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const SCALE_MAT<VECTOR<T,d> >& a,const SCALE_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(a.x-b.x);}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > operator- (const ID_MAT<VECTOR<T,d> >& a,const ID_MAT<VECTOR<T,d> >& b) {return SCALE_MAT<VECTOR<T,d> >(0);} // TODO: return ZERO_MAT

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

template<class T,int d> ZERO_MAT<VECTOR<T,d> > Times_Self_Transpose(const ZERO_MAT<VECTOR<T,d> >& a) {return a;}
template<class T,int d> ID_MAT<VECTOR<T,d> > Times_Self_Transpose(const ID_MAT<VECTOR<T,d> >& a) {return a;}
template<class T,int d> SCALE_MAT<VECTOR<T,d> > Times_Self_Transpose(const SCALE_MAT<VECTOR<T,d> >& a) {return a*a;}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Times_Self_Transpose(const MATRIX<T,d>& a) {return a.Outer_Product_Matrix();}
template<class T,int d> SYMMETRIC_MATRIX<T,d> Times_Self_Transpose(const SYMMETRIC_MATRIX<T,d>& a) {return a*a;}

template<int m,int s> M_FLAGS<m,s> Tst(M_FLAGS<m,s>);

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
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(ID_MAT<VECTOR<T,d> > z){return SYMMETRIC_MATRIX<T,d>()+1;}
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(SYMMETRIC_MATRIX<T,d> z){return z;}
template<class T,int d> inline MATRIX<T,d> Cast_Helper(const MATRIX<T,d>& z){return z;}

template<class TV,class H> struct CH_MAT;
template<class TV> struct CH_MAT<TV,M_FLAGS<0,0> > {typedef ZERO_MAT<TV> TYPE;};
template<class TV> struct CH_MAT<TV,M_FLAGS<1,0> > {typedef MATRIX<typename TV::SCALAR,TV::m> TYPE;};
template<class TV> struct CH_MAT<TV,M_FLAGS<0,1> > {typedef ID_MAT<TV> TYPE;};
template<class TV> struct CH_MAT<TV,M_FLAGS<1,1> > {typedef SCALE_MAT<TV> TYPE;};

template<class TV,class H> struct CH_SYM_MAT;
template<class TV> struct CH_SYM_MAT<TV,M_FLAGS<0,0> > {typedef ZERO_MAT<TV> TYPE;};
template<class TV> struct CH_SYM_MAT<TV,M_FLAGS<1,0> > {typedef SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> TYPE;};
template<class TV> struct CH_SYM_MAT<TV,M_FLAGS<0,1> > {typedef ID_MAT<TV> TYPE;};
template<class TV> struct CH_SYM_MAT<TV,M_FLAGS<1,1> > {typedef SCALE_MAT<TV> TYPE;};

#define MK_MAT(TV,H) typename CH_MAT<TV,decltype(H)>::TYPE
#define MK_SYM_MAT(TV,H) typename CH_SYM_MAT<TV,decltype(H)>::TYPE

template<class T,int d> void Fill_From_Helper(MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z){m=MATRIX<T,d>();}
template<class T,int d> void Fill_From_Helper(MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i){m=MATRIX<T,d>()+1;}
template<class T,int d> void Fill_From_Helper(MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s){m=MATRIX<T,d>()+s.x;}
template<class T,int d> void Fill_From_Helper(SYMMETRIC_MATRIX<T,d>& m,const ZERO_MAT<VECTOR<T,d> >& z){m=SYMMETRIC_MATRIX<T,d>();}
template<class T,int d> void Fill_From_Helper(SYMMETRIC_MATRIX<T,d>& m,const ID_MAT<VECTOR<T,d> >& i){m=SYMMETRIC_MATRIX<T,d>()+1;}
template<class T,int d> void Fill_From_Helper(SYMMETRIC_MATRIX<T,d>& m,const SCALE_MAT<VECTOR<T,d> >& s){m=SYMMETRIC_MATRIX<T,d>()+s.x;}
template<class T,int d> void Fill_From_Helper(SCALE_MAT<VECTOR<T,d> >& s,const ZERO_MAT<VECTOR<T,d> >& z) {s.x=0;}
template<class T,int d> void Fill_From_Helper(SCALE_MAT<VECTOR<T,d> >& s,const ID_MAT<VECTOR<T,d> >& i) {s.x=1;}


template<class T,int d> void
Outer_Product_Helper(MATRIX<T,d>& r,const VECTOR<T,d>& u,const VECTOR<T,d>& v)
{r=MATRIX<T,d>::Outer_Product(u,v);}

template<class TV,class TV2,class TV3> void
Outer_Product_Helper(ZERO_MAT<TV>& r,TV2 u,TV3 v) {}

template<class T,int d> void
Symmetric_Outer_Product_Helper(SYMMETRIC_MATRIX<T,d>& r,const VECTOR<T,d>& u,const VECTOR<T,d>& v)
{r=SYMMETRIC_MATRIX<T,d>::Symmetric_Outer_Product(u,v);}

template<class TV,class TV2,class TV3> void
Symmetric_Outer_Product_Helper(ZERO_MAT<TV>& r,TV2 u,TV3 v) {}

template<class T,int d> void
Outer_Product_Helper(SYMMETRIC_MATRIX<T,d>& r,const VECTOR<T,d>& u)
{r=SYMMETRIC_MATRIX<T,d>::Outer_Product(u);}

template<class TV,class TV2> void
Outer_Product_Helper(ZERO_MAT<TV>& r,TV2 u) {}
}
}
#endif
