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

struct DIFF_UNUSED;

template<class T,class VEC_IN>
struct GRADIENT
{
    static_assert(is_scalar<T>::value,"This definition of GRADIENT is only for scalars");
    static_assert(ASSERT_VALID_BLOCK_TYPES_VEC<VEC_IN>::value,"GRADIENT object is constructed from inconsistent block types");
    static_assert(!is_same<VEC_IN,DIFF_UNUSED>::value,"GRADIENT DIFF_UNUSED");
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
};

template<class T,class VEC,class VEC1>
void Fill_From(GRADIENT<T,VEC>& out,const GRADIENT<T,VEC1>& in)
{Fill_From(out.x,in.x);}

template<class T,class VEC>
auto operator* (T a,const GRADIENT<T,VEC>& u) -> decltype(u*a)
{return u*a;}

template<int i,class T,class A,int d,class VEC> inline void
Extract(VECTOR<A,d>& dx,const GRADIENT<T,VEC>& v)
{Extract<i>(dx,v.x);}

template<int i,class T,class OUT,class VEC>
void Get(OUT& o,const GRADIENT<T,VEC>& g)
{
    GET_VEC_HELPER<i>::f(o,g.x);
}

}
}
#endif
