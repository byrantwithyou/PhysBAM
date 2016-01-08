//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRADIENT_VEC
//##################################################################### 
#ifndef __GRADIENT_VEC__
#define __GRADIENT_VEC__

#include <Tools/Auto_Diff/GRADIENT.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV,class VEC>
struct GRADIENT_VEC
{
    typedef typename TV::SCALAR T;
    static_assert(ASSERT_VALID_BLOCK_TYPES_VEC<VEC>::value,"GRADIENT_VEC object is constructed from inconsistent block types");
    static_assert(!is_same<VEC,DIFF_UNUSED>::value,"GRADIENT_VEC DIFF_UNUSED");
    VEC x;

    GRADIENT_VEC operator+ () const {return *this;}
    auto operator- () const {GRADIENT_VEC<TV,decltype(VEC_NEG_B::Type(VEC()))> r;VEC_NEG_B()(r.x,x);return r;}
    auto operator* (T a) const {GRADIENT_VEC<TV,decltype(VEC_MUL_BS::Type(x,a))> r;VEC_MUL_BS()(r.x,x,a);return r;}
    auto operator/ (T a) const {GRADIENT_VEC<TV,decltype(VEC_DIV_BS::Type(x,a))> r;VEC_DIV_BS()(r.x,x,a);return r;}

    template<class VEC1>
    auto operator+ (const GRADIENT_VEC<TV,VEC1>& z) const {GRADIENT_VEC<TV,decltype(VEC_ADD_BB::Type(x,z.x))> r;VEC_ADD_BB()(r.x,x,z.x);return r;}

    template<class VEC1>
    auto operator- (const GRADIENT_VEC<TV,VEC1>& z) const {GRADIENT_VEC<TV,decltype(VEC_SUB_BB::Type(x,z.x))> r;VEC_SUB_BB()(r.x,x,z.x);return r;}

    template<class TV2>
    auto Transpose_Times(const TV2& w) const
        -> typename enable_if<IS_VECTOR<TV2>::value,GRADIENT<T,decltype(VEC_TRANSPOSE_TIMES::Type(this->x,w))> >::type
    {GRADIENT<T,decltype(VEC_TRANSPOSE_TIMES::Type(x,w))> r;VEC_TRANSPOSE_TIMES()(r.x,x,w);return r;}

    template<class T_MAT>
    auto Transpose_Times(const T_MAT& w) const
        -> typename enable_if<IS_MATRIX<T_MAT>::value,GRADIENT<T,decltype(VEC_TRANSPOSE_TIMES::Type(this->x,w))> >::type
    {GRADIENT_VEC<TV,decltype(VEC_TRANSPOSE_TIMES::Type(x,w))> r;VEC_TRANSPOSE_TIMES()(r.x,x,w);return r;}
};

template<class TV,class VEC,class VEC1> void
Fill_From(GRADIENT_VEC<TV,VEC>& out,const GRADIENT_VEC<TV,VEC1>& in)
{Fill_From(out.x,in.x);}

template<class TV,class VEC> auto
operator* (typename TV::SCALAR a,const GRADIENT_VEC<TV,VEC>& u){return u*a;}

template<class T,int d,class VEC> auto
Outer_Product(const VECTOR<T,d>& u,const GRADIENT<T,VEC>& v)
{GRADIENT_VEC<VECTOR<T,d>,decltype(VEC_OUTER_PRODUCT_REV::Type(VEC(),VECTOR<T,d>()))> r;VEC_OUTER_PRODUCT_REV()(r.x,v.x,u);return r;}

template<class TV,class VEC,class T_MAT> typename enable_if<IS_MATRIX<T_MAT>::value,GRADIENT_VEC<VECTOR<typename TV::SCALAR,T_MAT::m>,decltype(VEC_MUL_MB::Type(VEC(),T_MAT()))> >::type
operator* (const T_MAT& a,const GRADIENT_VEC<TV,VEC>& u)
{GRADIENT_VEC<VECTOR<typename TV::SCALAR,T_MAT::m>,decltype(VEC_MUL_MB::Type(VEC(),T_MAT()))> r;VEC_MUL_MB()(r.x,u.x,a);return r;}

template<int i,class TV,class A,int d,class VEC> inline void
Extract(VECTOR<A,d>& dx,const GRADIENT_VEC<TV,VEC>& v)
{Extract<i>(dx,v.x);}

template<int i,class TV,class OUT,class VEC>
void Get(OUT& o,const GRADIENT_VEC<TV,VEC>& g)
{
    GET_VEC_HELPER<i>::f(o,g.x);
}

}
}
#endif
