//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HESSIAN_VEC
//##################################################################### 
#ifndef __HESSIAN_VEC__
#define __HESSIAN_VEC__

#include <Tools/Auto_Diff/GRADIENT_VEC.h>
#include <Tools/Auto_Diff/HESSIAN.h>
#include <Tools/Auto_Diff/PRIMITIVE_MATRICES.h>
#include <Tools/Auto_Diff/PRIMITIVE_TENSORS.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
template<class TV,class MAT>
struct HESSIAN_VEC
{
    typedef typename TV::SCALAR T;
    MAT x;

    HESSIAN_VEC operator+ () const {return *this;}
    HESSIAN_VEC<TV,decltype(MAT_NEG::Type(MAT()))> operator- () const {HESSIAN_VEC<TV,decltype(MAT_NEG::Type(MAT()))> r;MAT_NEG()(r.x,x);return r;}
    HESSIAN_VEC<TV,decltype(MAT_SCALE::Type(MAT(),T()))> operator* (T a) const {HESSIAN_VEC<TV,decltype(MAT_SCALE::Type(MAT(),T()))> r;MAT_SCALE()(r.x,x,a);return r;}
    HESSIAN_VEC<TV,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> operator/ (T a) const {HESSIAN_VEC<TV,decltype(MAT_SCALE_DIV::Type(MAT(),T()))> r;MAT_SCALE_DIV()(r.x,x,a);return r;}

    template<class MAT2>
    HESSIAN_VEC<TV,decltype(MAT_ADD::Type(MAT(),MAT2()))> operator+ (const HESSIAN_VEC<TV,MAT2>& z) const {HESSIAN_VEC<TV,decltype(MAT_ADD::Type(MAT(),MAT2()))> r;MAT_ADD()(r.x,x,z.x);return r;}

    template<class MAT2>
    HESSIAN_VEC<TV,decltype(MAT_SUB::Type(MAT(),MAT2()))> operator- (const HESSIAN_VEC<TV,MAT2>& z) const {HESSIAN_VEC<TV,decltype(MAT_SUB::Type(MAT(),MAT2()))> r;MAT_SUB()(r.x,x,z.x);return r;}
};

template<class TV,class MAT,class MAT2>
void Fill_From(HESSIAN_VEC<TV,MAT>& out,const HESSIAN_VEC<TV,MAT2>& in)
{Fill_From(out.x,in.x);}

template<class TV,class MAT>
HESSIAN_VEC<TV,decltype(MAT_SCALE::Type(MAT(),typename TV::SCALAR()))> operator* (typename TV::SCALAR a,const HESSIAN_VEC<TV,MAT>& h){return h*a;}

template<class TV,class MAT,class TV2>
HESSIAN_VEC<TV,decltype(MAT_TENSOR_PRODUCT_0::Type(MAT(),TV2()))> Tensor_Product_0(const HESSIAN<TV,MAT>& m,const TV2& v)
{HESSIAN_VEC<TV,decltype(MAT_TENSOR_PRODUCT_0::Type(MAT(),TV2()))> r;MAT_TENSOR_PRODUCT_0()(r.x,m.x,v);return r;}

template<class TV,class MAT,class TV2>
HESSIAN_VEC<TV,decltype(MAT_TENSOR_PRODUCT_0::Type(MAT(),TV2()))> Symmetric_Tensor_Product_12(const HESSIAN<TV,MAT>& m,const TV2& v)
{HESSIAN_VEC<TV,decltype(MAT_TENSOR_PRODUCT_0::Type(MAT(),TV2()))> r;MAT_TENSOR_PRODUCT_0()(r.x,m.x,v);return r;}

template<class TV,class VEC,class VEC2>
HESSIAN_VEC<TV,decltype(MAT_SYM_TENSOR_PRODUCT_12::Type(VEC(),VEC2()))> Symmetric_Tensor_Product_12(const GRADIENT_VEC<TV,VEC>& gv,const GRADIENT<TV,VEC2>& g)
{HESSIAN_VEC<TV,decltype(MAT_SYM_TENSOR_PRODUCT_12::Type(VEC(),VEC2()))> r;MAT_SYM_TENSOR_PRODUCT_12()(r.x,gv.x,g.x);return r;}

template<class TV,class MAT,class TV2>
typename ENABLE_IF<IS_VECTOR<TV2>::value,HESSIAN<TV,decltype(MAT_CONTRACT_0::Type(MAT(),TV2()))> >::TYPE Contract_0(const HESSIAN_VEC<TV,MAT>& h,const TV2& v)
{HESSIAN<TV,decltype(MAT_CONTRACT_0::Type(MAT(),TV2()))> r;MAT_CONTRACT_0()(r.x,h.x,v);return r;}

template<class TV,class MAT,class T_MAT>
typename ENABLE_IF<IS_MATRIX<T_MAT>::value,HESSIAN_VEC<TV,decltype(MAT_CONTRACT_0::Type(MAT(),T_MAT()))> >::TYPE Contract_0(const HESSIAN_VEC<TV,MAT>& h,const T_MAT& v)
{HESSIAN_VEC<TV,decltype(MAT_CONTRACT_0::Type(MAT(),T_MAT()))> r;MAT_CONTRACT_0()(r.x,h.x,v);return r;}

template<class TV,class VEC,class VEC2>
HESSIAN_VEC<TV,decltype(MAT_SYM_DOUBLE_CONTRACT_12::Type(VEC(),VEC2(),PERM_TENSOR<TV>()))> 
Symmetric_Double_Contract_12_With_Tensor(const PERM_TENSOR<TV>& t,const GRADIENT_VEC<TV,VEC>& a,const GRADIENT_VEC<TV,VEC2>& b)
{HESSIAN_VEC<TV,decltype(MAT_SYM_DOUBLE_CONTRACT_12::Type(VEC(),VEC2(),PERM_TENSOR<TV>()))> r;MAT_SYM_DOUBLE_CONTRACT_12()(r.x,a.x,b.x,t);return r;}
}
}


#endif
