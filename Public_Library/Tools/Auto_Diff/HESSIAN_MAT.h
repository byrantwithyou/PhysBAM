//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HESSIAN_MAT
//##################################################################### 
#ifndef __HESSIAN_MAT__
#define __HESSIAN_MAT__

#include <Tools/Auto_Diff/GRADIENT_MAT.h>
#include <Tools/Auto_Diff/HESSIAN_VEC.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
template<class T_MAT,class MAT>
struct HESSIAN_MAT
{
    typedef typename T_MAT::SCALAR T;
    static_assert(!is_same<MAT,DIFF_UNUSED>::value,"HESSIAN_MAT DIFF_UNUSED");
    static_assert(ASSERT_VALID_BLOCK_TYPES_MAT<MAT>::value,"HESSIAN_MAT object is constructed from inconsistent block types");
    MAT x;

    HESSIAN_MAT operator+ () const {return *this;}
    HESSIAN_MAT<T_MAT,decltype(MAT_NEG::Type(MAT()))> operator- () const {HESSIAN_MAT<T_MAT,decltype(MAT_NEG::Type(MAT()))> r;MAT_NEG()(r.x,x);return r;}
    HESSIAN_MAT<T_MAT,decltype(MAT_SCALE::Type(MAT(),T()))> operator* (T a) const {HESSIAN_MAT<T_MAT,decltype(MAT_SCALE::Type(MAT(),T()))> r;MAT_SCALE()(r.x,x,a);return r;}
    HESSIAN_MAT<T_MAT,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> operator/ (T a) const {HESSIAN_MAT<T_MAT,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> r;MAT_SCALE_DIV()(r.x,x,a);return r;}

    template<class MAT1,class T_MAT1>
    HESSIAN_MAT<T_MAT,decltype(MAT_ADD::Type(MAT(),MAT1()))> operator+ (const HESSIAN_MAT<T_MAT1,MAT1>& z) const {HESSIAN_MAT<decltype(x+z.x),decltype(MAT_ADD::Type(MAT(),MAT1()))> r;MAT_ADD()(r.x,x,z.x);return r;}

    template<class MAT1,class T_MAT1>
    HESSIAN_MAT<T_MAT,decltype(MAT_SUB::Type(MAT(),MAT1()))> operator- (const HESSIAN_MAT<T_MAT1,MAT1>& z) const {HESSIAN_MAT<decltype(x-z.x),decltype(MAT_SUB::Type(MAT(),MAT1()))> r;MAT_SUB()(r.x,x,z.x);return r;}
};

template<class T_MAT,class T_MAT1,class MAT,class MAT1>
void Fill_From(HESSIAN_MAT<T_MAT,MAT>& out,const HESSIAN_MAT<T_MAT1,MAT1>& in)
{Fill_From(out.x,in.x);}

template<int i,int j,class T_MAT,class A,int d,class MAT> inline void
Extract(MATRIX<A,d>& ddx,const HESSIAN_MAT<T_MAT,MAT>& h)
{Extract<i,j>(ddx,h.x);}

template<int i,int j,class T_MAT,class OUT,class MAT>
void Get(OUT& o,const HESSIAN_MAT<T_MAT,MAT>& h)
{
    GET_MAT_HELPER<i,j>::f(o,h.x);
}

template<class T_MAT,class MAT>
auto operator* (typename T_MAT::SCALAR a,const HESSIAN_MAT<T_MAT,MAT>& h) -> decltype(h*a) {return h*a;}
}
}
#endif
