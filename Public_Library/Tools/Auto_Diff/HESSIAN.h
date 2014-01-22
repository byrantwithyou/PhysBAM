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
#include <Tools/Auto_Diff/PRIMITIVE_MATRICES.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV,class MAT>
struct HESSIAN
{
    typedef typename TV::SCALAR T;
    MAT x;

    HESSIAN operator+ () const {return *this;}
    HESSIAN<TV,decltype(MAT_NEG::Type(MAT()))> operator- () const {HESSIAN<TV,decltype(MAT_NEG::Type(MAT()))> r;MAT_NEG()(r.x,x);return r;}
    HESSIAN<TV,decltype(MAT_SCALE::Type(MAT(),T()))> operator* (T a) const {HESSIAN<TV,decltype(MAT_SCALE::Type(MAT(),T()))> r;MAT_SCALE()(r.x,x,a);return r;}
    HESSIAN<TV,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> operator/ (T a) const {HESSIAN<TV,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> r;MAT_SCALE_DIV()(r.x,x,a);return r;}

    template<class MAT1>
    HESSIAN<TV,decltype(MAT_ADD::Type(MAT(),MAT1()))> operator+ (const HESSIAN<TV,MAT1>& z) const {HESSIAN<TV,decltype(MAT_ADD::Type(MAT(),MAT1()))> r;MAT_ADD()(r.x,x,z.x);return r;}

    template<class MAT1>
    HESSIAN<TV,decltype(MAT_SUB::Type(MAT(),MAT1()))> operator- (const HESSIAN<TV,MAT1>& z) const {HESSIAN<TV,decltype(MAT_SUB::Type(MAT(),MAT1()))> r;MAT_SUB()(r.x,x,z.x);return r;}

    MATRIX<T,TV::m> operator()(int i,int j) const {MATRIX<T,TV::m> m;Get(m,x,i,j);return m;}

    SYMMETRIC_MATRIX<T,TV::m> operator()(int i) const {SYMMETRIC_MATRIX<T,TV::m> m;Get_Diag(m,x,i);return m;}
};

template<class TV,class MAT,class MAT1>
void Fill_From(HESSIAN<TV,MAT>& out,const HESSIAN<TV,MAT1>& in)
{Fill_From(out.x,in.x);}

template<class TV,class MAT>
HESSIAN<TV,decltype(MAT_SCALE::Type(MAT(),typename TV::SCALAR()))> operator* (typename TV::SCALAR a,const HESSIAN<TV,MAT>& h){return h*a;}

template<class TV,class VEC,class VEC1> HESSIAN<TV,decltype(MAT_SYM_OUTER_PRODUCT::Type(VEC(),VEC1()))>
Symmetric_Outer_Product(const GRADIENT<TV,VEC>& u,const GRADIENT<TV,VEC1>& v)
{HESSIAN<TV,decltype(MAT_SYM_OUTER_PRODUCT::Type(VEC(),VEC1()))> H;MAT_SYM_OUTER_PRODUCT()(H.x,u.x,v.x);return H;}

template<class TV,class VEC> HESSIAN<TV,decltype(MAT_OUTER_PRODUCT::Type(VEC()))>
Outer_Product(const GRADIENT<TV,VEC>& u)
{HESSIAN<TV,decltype(MAT_OUTER_PRODUCT::Type(VEC()))> H;MAT_OUTER_PRODUCT()(H.x,u.x);return H;}

template<class TV,class VEC,class VEC1> HESSIAN<TV,decltype(MAT_SYM_TRANSPOSE_TIMES::Type(VEC(),VEC1()))>
Symmetric_Transpose_Times(const GRADIENT_VEC<TV,VEC>& u,const GRADIENT_VEC<TV,VEC1>& v)
{HESSIAN<TV,decltype(MAT_SYM_TRANSPOSE_TIMES::Type(VEC(),VEC1()))> H;MAT_SYM_TRANSPOSE_TIMES()(H.x,u.x,v.x);return H;}

template<class TV,class VEC> HESSIAN<TV,decltype(MAT_TRANSPOSE_TIMES_SELF::Type(VEC()))>
Transpose_Times_Self(const GRADIENT_VEC<TV,VEC>& u)
{HESSIAN<TV,decltype(MAT_TRANSPOSE_TIMES_SELF::Type(VEC()))> H;MAT_TRANSPOSE_TIMES_SELF()(H.x,u.x);return H;}

}
}


#endif
