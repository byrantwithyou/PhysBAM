//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRADIENT
//##################################################################### 
#ifndef __GRADIENT__
#define __GRADIENT__

#include <Tools/Auto_Diff/VEC_HOLDER.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class T,int m> struct GRADIENT_GET_HELPER {typedef VECTOR<T,m> V;};
template<class T> struct GRADIENT_GET_HELPER<T,-1> {typedef T V;};

template<class T,class VEC_IN>
struct GRADIENT
{
    static_assert(is_scalar<T>::value,"This definition of GRADIENT is only for scalars");
    typedef VEC_IN VEC;
    VEC x;

    GRADIENT operator+ () const {return *this;}
    GRADIENT<T,decltype(VEC_NEG::Type(VEC()))> operator- () const {GRADIENT<T,decltype(VEC_NEG::Type(VEC()))> r;VEC_NEG()(r.x,x);return r;}
    GRADIENT<T,decltype(VEC_SCALE::Type(VEC(),T()))> operator* (T a) const {GRADIENT<T,decltype(VEC_SCALE::Type(VEC(),T()))> r;VEC_SCALE()(r.x,x,a);return r;}
    GRADIENT<T,decltype(VEC_SCALE_DIV::Type(VEC(),T()))> operator/ (T a) const {GRADIENT<T,decltype(VEC_SCALE_DIV::Type(VEC(),T()))> r;VEC_SCALE_DIV()(r.x,x,a);return r;}

    template<class VEC1>
    GRADIENT<T,decltype(VEC_ADD::Type(VEC(),VEC1()))> operator+ (const GRADIENT<T,VEC1>& z) const {GRADIENT<T,decltype(VEC_ADD::Type(VEC(),VEC1()))> r;VEC_ADD()(r.x,x,z.x);return r;}

    template<class VEC1>
    GRADIENT<T,decltype(VEC_SUB::Type(VEC(),VEC1()))> operator- (const GRADIENT<T,VEC1>& z) const {GRADIENT<T,decltype(VEC_SUB::Type(VEC(),VEC1()))> r;VEC_SUB()(r.x,x,z.x);return r;}

    template<int m>
    typename GRADIENT_GET_HELPER<T,m>::V Get_Block(int i) const {typename GRADIENT_GET_HELPER<T,m>::V v;Get(v,x,i);return v;}

    template<int i,int m>
    void Set_Block(const VECTOR<T,m>& v) const {Set<i>(v,x);}

    template<int i>
    void Set_Block(const T& v) const {Set<i>(v,x);}
};

template<class T,class VEC,class VEC1>
void Fill_From(GRADIENT<T,VEC>& out,const GRADIENT<T,VEC1>& in)
{Fill_From(out.x,in.x);}

template<class T,class VEC>
auto operator* (T a,const GRADIENT<T,VEC>& u) -> decltype(u*a)
{return u*a;}

template<class T,int m,int d,class VEC> inline void
Extract(VECTOR<VECTOR<T,m>,d>& dx,const GRADIENT<T,VEC>& v)
{Extract(dx,v.x);}

template<class T,int m,int d,class VEC> inline void
Extract(MATRIX<T,m,d>& dx,const GRADIENT<T,VEC>& v)
{VECTOR<T,m> w;for(int i=0;i<d;i++){Get(w,v.x,i);dx.Set_Column(i,w);}}

}
}
#endif
