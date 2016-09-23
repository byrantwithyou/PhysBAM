//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HESSIAN_MAT
//##################################################################### 
#ifndef __HESSIAN_MAT__
#define __HESSIAN_MAT__

#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/GRADIENT_MAT.h>
#include <Tools/Auto_Diff/HESSIAN_VEC.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
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

    auto operator+ () const {return *this;}
    auto operator- () const {HESSIAN_MAT<T_MAT,decltype(MAT_NEG_B::Type(x))> r;MAT_NEG_B()(r.x,x);return r;}
    auto operator* (T a) const {HESSIAN_MAT<T_MAT,decltype(MAT_MUL_BS::Type(x,a))> r;MAT_MUL_BS()(r.x,x,a);return r;}
    auto operator/ (T a) const {HESSIAN_MAT<T_MAT,decltype(MAT_DIV_BS::Type(x,a))> r;MAT_DIV_BS()(r.x,x,a);return r;}

    template<class MAT1,class T_MAT1>
    auto operator+ (const HESSIAN_MAT<T_MAT1,MAT1>& z) const {HESSIAN_MAT<decltype(x+z.x),decltype(MAT_ADD_BB::Type(x,z.x))> r;MAT_ADD_BB()(r.x,x,z.x);return r;}

    template<class MAT1,class T_MAT1>
    auto operator- (const HESSIAN_MAT<T_MAT1,MAT1>& z) const {HESSIAN_MAT<decltype(x-z.x),decltype(MAT_SUB_BB::Type(x,z.x))> r;MAT_SUB_BB()(r.x,x,z.x);return r;}
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
auto operator* (typename T_MAT::SCALAR a,const HESSIAN_MAT<T_MAT,MAT>& h) {return h*a;}
}
}
#endif
