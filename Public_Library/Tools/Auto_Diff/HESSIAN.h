//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HESSIAN
//##################################################################### 
#ifndef __HESSIAN__
#define __HESSIAN__

#include <Tools/Auto_Diff/GRADIENT.h>
#include <Tools/Auto_Diff/GRADIENT_VEC.h>
#include <Tools/Auto_Diff/MAT_HOLDER.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class T,int m,int n> struct HESSIAN_GET_HELPER {typedef MATRIX<T,m,n> M;typedef SYMMETRIC_MATRIX<T,m> SM;};
template<class T,int m> struct HESSIAN_GET_HELPER<T,-1,m> {typedef VECTOR<T,m> M;};
template<class T,int m> struct HESSIAN_GET_HELPER<T,m,-1> {typedef VECTOR<T,m> M;};
template<class T> struct HESSIAN_GET_HELPER<T,-1,-1> {typedef T M;typedef T SM;};

template<class T,class MAT>
struct HESSIAN
{
    static_assert(is_scalar<T>::value,"This definition of HESSIAN is only for scalars");
    static_assert(ASSERT_VALID_BLOCK_TYPES_MAT<MAT>::value,"HESSIAN object is constructed from inconsistent block types");
    static_assert(!is_same<MAT,DIFF_UNUSED>::value,"HESSIAN DIFF_UNUSED");
    MAT x;

    HESSIAN operator+ () const {return *this;}
    HESSIAN<T,decltype(MAT_NEG::Type(MAT()))> operator- () const {HESSIAN<T,decltype(MAT_NEG::Type(MAT()))> r;MAT_NEG()(r.x,x);return r;}
    HESSIAN<T,decltype(MAT_SCALE::Type(MAT(),T()))> operator* (T a) const {HESSIAN<T,decltype(MAT_SCALE::Type(MAT(),T()))> r;MAT_SCALE()(r.x,x,a);return r;}
    HESSIAN<T,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> operator/ (T a) const {HESSIAN<T,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> r;MAT_SCALE_DIV()(r.x,x,a);return r;}

    template<class MAT1>
    HESSIAN<T,decltype(MAT_ADD::Type(MAT(),MAT1()))> operator+ (const HESSIAN<T,MAT1>& z) const {HESSIAN<T,decltype(MAT_ADD::Type(MAT(),MAT1()))> r;MAT_ADD()(r.x,x,z.x);return r;}

    template<class MAT1>
    HESSIAN<T,decltype(MAT_SUB::Type(MAT(),MAT1()))> operator- (const HESSIAN<T,MAT1>& z) const {HESSIAN<T,decltype(MAT_SUB::Type(MAT(),MAT1()))> r;MAT_SUB()(r.x,x,z.x);return r;}
};

template<class T,class MAT,class MAT1>
void Fill_From(HESSIAN<T,MAT>& out,const HESSIAN<T,MAT1>& in)
{Fill_From(out.x,in.x);}

template<int i,int j,class T,class A,int d,class MAT> inline void
Extract(MATRIX<A,d>& ddx,const HESSIAN<T,MAT>& h)
{Extract<i,j>(ddx,h.x);}

template<int i,int j,class T,class OUT,class MAT>
void Get(OUT& o,const HESSIAN<T,MAT>& h)
{
    GET_MAT_HELPER<i,j>::f(o,h.x);
}

template<class T,class MAT>
HESSIAN<T,decltype(MAT_SCALE::Type(MAT(),T()))> operator* (T a,const HESSIAN<T,MAT>& h){return h*a;}

template<class T,class VEC,class VEC1> HESSIAN<T,decltype(MAT_SYM_OUTER_PRODUCT::Type(VEC(),VEC1()))>
Symmetric_Outer_Product(const GRADIENT<T,VEC>& u,const GRADIENT<T,VEC1>& v)
{HESSIAN<T,decltype(MAT_SYM_OUTER_PRODUCT::Type(VEC(),VEC1()))> H;MAT_SYM_OUTER_PRODUCT()(H.x,u.x,v.x);return H;}

template<class T,class VEC> HESSIAN<T,decltype(MAT_OUTER_PRODUCT::Type(VEC()))>
Outer_Product(const GRADIENT<T,VEC>& u)
{HESSIAN<T,decltype(MAT_OUTER_PRODUCT::Type(VEC()))> H;MAT_OUTER_PRODUCT()(H.x,u.x);return H;}

template<class TV,class VEC,class VEC1> HESSIAN<typename TV::SCALAR,decltype(MAT_SYM_TRANSPOSE_TIMES::Type(VEC(),VEC1()))>
Symmetric_Transpose_Times(const GRADIENT_VEC<TV,VEC>& u,const GRADIENT_VEC<TV,VEC1>& v)
{HESSIAN<typename TV::SCALAR,decltype(MAT_SYM_TRANSPOSE_TIMES::Type(VEC(),VEC1()))> H;MAT_SYM_TRANSPOSE_TIMES()(H.x,u.x,v.x);return H;}

template<class TV,class VEC> HESSIAN<typename TV::SCALAR,decltype(MAT_TRANSPOSE_TIMES_SELF::Type(VEC()))>
Transpose_Times_Self(const GRADIENT_VEC<TV,VEC>& u)
{HESSIAN<typename TV::SCALAR,decltype(MAT_TRANSPOSE_TIMES_SELF::Type(VEC()))> H;MAT_TRANSPOSE_TIMES_SELF()(H.x,u.x);return H;}

}
}


#endif
