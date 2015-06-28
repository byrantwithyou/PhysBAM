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
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
template<class TV,class MAT>
struct HESSIAN_VEC
{
    typedef typename TV::SCALAR T;
    static_assert(!is_same<MAT,DIFF_UNUSED>::value,"HESSIAN_VEC DIFF_UNUSED");
    static_assert(ASSERT_VALID_BLOCK_TYPES_MAT<MAT>::value,"HESSIAN_VEC object is constructed from inconsistent block types");
    MAT x;

    auto operator+ () const {return *this;}
    auto operator- () const {HESSIAN_VEC<TV,decltype(MAT_NEG_B::Type(x))> r;MAT_NEG_B()(r.x,x);return r;}
    auto operator* (T a) const {HESSIAN_VEC<TV,decltype(MAT_MUL_BS::Type(x,T()))> r;MAT_MUL_BS()(r.x,x,a);return r;}
    auto operator/ (T a) const {HESSIAN_VEC<TV,decltype(MAT_DIV_BS::Type(x,T()))> r;MAT_DIV_BS()(r.x,x,a);return r;}

    template<class MAT1>
    auto operator+ (const HESSIAN_VEC<TV,MAT1>& z) const {HESSIAN_VEC<TV,decltype(MAT_ADD_BB::Type(x,z.x))> r;MAT_ADD_BB()(r.x,x,z.x);return r;}

    template<class MAT1>
    auto operator- (const HESSIAN_VEC<TV,MAT1>& z) const {HESSIAN_VEC<TV,decltype(MAT_SUB_BB::Type(x,z.x))> r;MAT_SUB_BB()(r.x,x,z.x);return r;}
};

template<class TV,class MAT,class MAT1>
void Fill_From(HESSIAN_VEC<TV,MAT>& out,const HESSIAN_VEC<TV,MAT1>& in)
{Fill_From(out.x,in.x);}

template<int i,int j,class TV,class A,int d,class MAT> inline void
Extract(MATRIX<A,d>& ddx,const HESSIAN_VEC<TV,MAT>& h)
{Extract<i,j>(ddx,h.x);}

template<int i,int j,class TV,class OUT,class MAT>
void Get(OUT& o,const HESSIAN_VEC<TV,MAT>& h)
{
    GET_MAT_HELPER<i,j>::f(o,h.x);
}

template<class TV,class MAT> auto
operator* (typename TV::SCALAR a,const HESSIAN_VEC<TV,MAT>& h){return h*a;}

template<class TV2,class MAT> auto
Tensor_Product_0(const HESSIAN<typename TV2::SCALAR,MAT>& m,const TV2& v)
{HESSIAN_VEC<VECTOR<typename TV2::SCALAR,TV2::m>,decltype(MAT_TENSOR_PRODUCT_0::Type(MAT(),TV2()))> r;MAT_TENSOR_PRODUCT_0()(r.x,m.x,v);return r;}

template<class MAT,class TV2> auto
Symmetric_Tensor_Product_12(const HESSIAN<typename TV2::SCALAR,MAT>& m,const TV2& v)
{HESSIAN_VEC<VECTOR<typename TV2::SCALAR,TV2::m>,decltype(MAT_TENSOR_PRODUCT_0::Type(MAT(),TV2()))> r;MAT_TENSOR_PRODUCT_0()(r.x,m.x,v);return r;}

template<class T,class TV,class VEC,class VEC1> auto
Symmetric_Tensor_Product_12(const GRADIENT_VEC<TV,VEC>& gv,const GRADIENT<T,VEC1>& g)
{HESSIAN_VEC<TV,decltype(MAT_SYM_TENSOR_PRODUCT_12::Type(VEC(),VEC1()))> r;MAT_SYM_TENSOR_PRODUCT_12()(r.x,gv.x,g.x);return r;}

template<class TV,class MAT,class TV2> auto
Contract_0(const HESSIAN_VEC<TV,MAT>& h,const TV2& v,typename enable_if<IS_VECTOR<TV2>::value,void*>::type=0)
{HESSIAN<typename TV::SCALAR,decltype(MAT_CONTRACT_0::Type(MAT(),TV2()))> r;MAT_CONTRACT_0()(r.x,h.x,v);return r;}

template<class TV,class MAT,class T_MAT> auto
Contract_00(const HESSIAN_VEC<TV,MAT>& h,const T_MAT& v,typename enable_if<IS_MATRIX<T_MAT>::value,void*>::type=0)
{HESSIAN_VEC<TV,decltype(MAT_CONTRACT_00::Type(MAT(),T_MAT()))> r;MAT_CONTRACT_00()(r.x,h.x,v);return r;}

template<class TV,class VEC,class VEC1,class T_TEN> auto
Symmetric_Double_Contract_12_With_Tensor(const T_TEN& t,const GRADIENT_VEC<TV,VEC>& a,const GRADIENT_VEC<TV,VEC1>& b)
{HESSIAN_VEC<TV,decltype(MAT_SYM_DOUBLE_CONTRACT_12::Type(VEC(),VEC1(),T_TEN()))> r;MAT_SYM_DOUBLE_CONTRACT_12()(r.x,a.x,b.x,t);return r;}

template<class TV,class MAT,class T_MAT> auto
Contract_01(const HESSIAN_VEC<TV,MAT>& h,const T_MAT& v,typename enable_if<IS_MATRIX<T_MAT>::value,void*>::type=0)
{HESSIAN_VEC<TV,decltype(MAT_CONTRACT_01::Type(MAT(),T_MAT()))> r;MAT_CONTRACT_01()(r.x,h.x,v);return r;}
}
}


#endif
