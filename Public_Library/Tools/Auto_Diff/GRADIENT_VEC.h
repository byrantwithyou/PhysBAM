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
    VEC x;

    GRADIENT_VEC operator+ () const {return *this;}
    GRADIENT_VEC<TV,decltype(VEC_NEG::Type(VEC()))> operator- () const {GRADIENT_VEC<TV,decltype(VEC_NEG::Type(VEC()))> r;VEC_NEG()(r.x,x);return r;}
    GRADIENT_VEC<TV,decltype(VEC_SCALE::Type(VEC(),T()))> operator* (T a) const {GRADIENT_VEC<TV,decltype(VEC_SCALE::Type(VEC(),T()))> r;VEC_SCALE()(r.x,x,a);return r;}
    GRADIENT_VEC<TV,decltype(VEC_SCALE_DIV::Type(VEC(),T()))> operator/ (T a) const {GRADIENT_VEC<TV,decltype(VEC_SCALE_DIV::Type(VEC(),T()))> r;VEC_SCALE_DIV()(r.x,x,a);return r;}

    template<class VEC1>
    GRADIENT_VEC<TV,decltype(VEC_ADD::Type(VEC(),VEC1()))> operator+ (const GRADIENT_VEC<TV,VEC1>& z) const {GRADIENT_VEC<TV,decltype(VEC_ADD::Type(VEC(),VEC1()))> r;VEC_ADD()(r.x,x,z.x);return r;}

    template<class VEC1>
    GRADIENT_VEC<TV,decltype(VEC_SUB::Type(VEC(),VEC1()))> operator- (const GRADIENT_VEC<TV,VEC1>& z) const {GRADIENT_VEC<TV,decltype(VEC_SUB::Type(VEC(),VEC1()))> r;VEC_SUB()(r.x,x,z.x);return r;}

    template<class TV2>
    typename enable_if<IS_VECTOR<TV2>::value,GRADIENT<T,decltype(VEC_TRANSPOSE_TIMES::Type(VEC(),TV2()))> >::type Transpose_Times(const TV2& w) const
    {GRADIENT<T,decltype(VEC_TRANSPOSE_TIMES::Type(VEC(),TV2()))> r;VEC_TRANSPOSE_TIMES()(r.x,x,w);return r;}

    template<class T_MAT>
    typename enable_if<IS_MATRIX<T_MAT>::value,GRADIENT_VEC<TV,decltype(VEC_TRANSPOSE_TIMES::Type(VEC(),T_MAT()))> >::type Transpose_Times(const T_MAT& w) const
    {GRADIENT_VEC<TV,decltype(VEC_TRANSPOSE_TIMES::Type(VEC(),T_MAT()))> r;VEC_TRANSPOSE_TIMES()(r.x,x,w);return r;}
};

template<class TV,class VEC,class VEC1>
void Fill_From(GRADIENT_VEC<TV,VEC>& out,const GRADIENT_VEC<TV,VEC1>& in)
{Fill_From(out.x,in.x);}

template<class TV,class VEC>
GRADIENT_VEC<TV,decltype(VEC_SCALE::Type(VEC(),typename TV::SCALAR()))> operator* (typename TV::SCALAR a,const GRADIENT_VEC<TV,VEC>& u){return u*a;}

template<class T,int d,class VEC> GRADIENT_VEC<VECTOR<T,d>,decltype(VEC_OUTER_PRODUCT_REV::Type(VEC(),VECTOR<T,d>()))>
Outer_Product(const VECTOR<T,d>& u,const GRADIENT<T,VEC>& v)
{GRADIENT_VEC<VECTOR<T,d>,decltype(VEC_OUTER_PRODUCT_REV::Type(VEC(),VECTOR<T,d>()))> r;VEC_OUTER_PRODUCT_REV()(r.x,v.x,u);return r;}

template<class TV,class VEC,class T_MAT> typename enable_if<IS_MATRIX<T_MAT>::value,GRADIENT_VEC<TV,decltype(VEC_SCALE_REV::Type(VEC(),T_MAT()))> >::type
operator* (const T_MAT& a,const GRADIENT_VEC<TV,VEC>& u)
{GRADIENT_VEC<TV,decltype(VEC_SCALE_REV::Type(VEC(),T_MAT()))> r;VEC_SCALE_REV()(r.x,u.x,a);return r;}
}
}
#endif
