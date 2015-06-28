//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRADIENT_MAT
//##################################################################### 
#ifndef __GRADIENT_MAT__
#define __GRADIENT_MAT__

#include <Tools/Auto_Diff/GRADIENT_VEC.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class T_MAT,class VEC>
struct GRADIENT_MAT
{
    typedef typename T_MAT::SCALAR T;
    static_assert(ASSERT_VALID_BLOCK_TYPES_VEC<VEC>::value,"GRADIENT_MAT object is constructed from inconsistent block types");
    static_assert(!is_same<VEC,DIFF_UNUSED>::value,"GRADIENT_MAT DIFF_UNUSED");
    static_assert(is_same<T_MAT,MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n> >::value,"make T_MAT be MATRIX<T,m,n>");
    VEC x;

    GRADIENT_MAT operator+ () const {return *this;}
    auto operator- () const {GRADIENT_MAT<T_MAT,decltype(VEC_NEG_B::Type(VEC()))> r;VEC_NEG_B()(r.x,x);return r;}
    auto operator* (T a) const {GRADIENT_MAT<T_MAT,decltype(VEC_MUL_BS::Type(VEC(),T()))> r;VEC_MUL_BS()(r.x,x,a);return r;}
    auto operator/ (T a) const {GRADIENT_MAT<T_MAT,decltype(VEC_DIV_BS::Type(VEC(),T()))> r;VEC_DIV_BS()(r.x,x,a);return r;}

    template<class VEC1,class T_MAT2> auto
    operator+ (const GRADIENT_MAT<T_MAT2,VEC1>& z) const {GRADIENT_MAT<T_MAT,decltype(VEC_ADD_BB::Type(x,z.x))> r;VEC_ADD_BB()(r.x,x,z.x);return r;}

    template<class VEC1,class T_MAT2> auto
    operator- (const GRADIENT_MAT<T_MAT2,VEC1>& z) const {GRADIENT_MAT<T_MAT,decltype(VEC_SUB_BB::Type(x,z.x))> r;VEC_SUB_BB()(r.x,x,z.x);return r;}
};

template<class T_MAT,class VEC,class T_MAT2,class VEC1>
void Fill_From(GRADIENT_MAT<T_MAT,VEC>& out,const GRADIENT_MAT<T_MAT2,VEC1>& in)
{Fill_From(out.x,in.x);}

template<class T_MAT,class VEC>
auto operator* (typename T_MAT::SCALAR a,const GRADIENT_MAT<T_MAT,VEC>& u) {return u*a;}

template<class TV,class T,int d,class VEC> auto
Tensor_Product_1(const GRADIENT_VEC<TV,VEC>& v,const VECTOR<T,d>& u)
{GRADIENT_MAT<MATRIX<T,TV::m,d>,decltype(TENSOR_PRODUCT_1_VM_V::Type(VEC(),u))> r;TENSOR_PRODUCT_1_VM_V()(r.x,v.x,u);return r;}

template<class TV,class T,int d,class VEC> auto
Tensor_Product_0(const GRADIENT_VEC<TV,VEC>& v,const VECTOR<T,d>& u)
{GRADIENT_MAT<MATRIX<T,TV::m,d>,decltype(TENSOR_PRODUCT_0_VM_V::Type(VEC(),u))> r;TENSOR_PRODUCT_0_VM_V()(r.x,v.x,u);return r;}

template<class T_MAT,class T,class VEC> auto
Tensor_Product_2(const T_MAT& m,const GRADIENT<T,VEC>& v)
    -> typename enable_if<IS_MATRIX<T_MAT>::value,GRADIENT_MAT<MATRIX<T,T_MAT::m,T_MAT::n>,decltype(TENSOR_PRODUCT_2_VV_M::Type(v.x,m))> >::type
{GRADIENT_MAT<MATRIX<T,T_MAT::m,T_MAT::n>,decltype(TENSOR_PRODUCT_2_VV_M::Type(v.x,m))> r;TENSOR_PRODUCT_2_VV_M()(r.x,v.x,m);return r;}

template<class T_MAT,class MAT> auto
Twice_Symmetric_Part_01(const GRADIENT_MAT<T_MAT,MAT>& v)
{GRADIENT_MAT<T_MAT,decltype(TWICE_SYMMETRIC_PART_01_VT::Type(v.x))> r;TWICE_SYMMETRIC_PART_01_VT()(r.x,v.x);return r;}

template<class T_MAT,class MAT> auto
Contract_01(const GRADIENT_MAT<T_MAT,MAT>& v)
{GRADIENT<typename T_MAT::SCALAR,decltype(CONTRACT_01_VT::Type(v.x))> r;CONTRACT_01_VT()(r.x,v.x);return r;}

template<class T_MAT,class MAT> auto
Transposed_01(const GRADIENT_MAT<T_MAT,MAT>& v)
{GRADIENT_MAT<decltype(T_MAT().Transposed()),decltype(TRANSPOSED_01_VT::Type(v.x))> r;TRANSPOSED_01_VT()(r.x,v.x);return r;}

template<class T_MAT,class MAT,class T_MAT2> auto
Contract_00(const GRADIENT_MAT<T_MAT,MAT>& v,const T_MAT2& m)
{GRADIENT_MAT<decltype(Transpose_Times(T_MAT(),T_MAT2())),decltype(CONTRACT_00_VT_M::Type(v.x,m))> r;CONTRACT_00_VT_M()(r.x,v.x,m);return r;}

template<class T_MAT,class MAT,class T_MAT2> auto
Contract_10(const GRADIENT_MAT<T_MAT,MAT>& v,const T_MAT2& m)
{GRADIENT_MAT<decltype(Transpose_Times(T_MAT(),T_MAT2())),decltype(CONTRACT_10_VT_M::Type(v.x,m))> r;CONTRACT_00_VT_M()(r.x,v.x,m);return r;}

template<class T_MAT,class MAT,class T_MAT2> auto
Contract_01(const GRADIENT_MAT<T_MAT,MAT>& v,const T_MAT2& m)
{GRADIENT_MAT<decltype(Transpose_Times(T_MAT(),T_MAT2())),decltype(CONTRACT_01_VT_M::Type(v.x,m))> r;CONTRACT_00_VT_M()(r.x,v.x,m);return r;}

template<class T_MAT,class MAT,class T_MAT2> auto
Contract_11(const GRADIENT_MAT<T_MAT,MAT>& v,const T_MAT2& m)
{GRADIENT_MAT<decltype(Transpose_Times(T_MAT(),T_MAT2())),decltype(CONTRACT_11_VT_M::Type(v.x,m))> r;CONTRACT_11_VT_M()(r.x,v.x,m);return r;}

template<class T_MAT,class MAT,class T_MAT2> auto
Double_Contract_00_11(const GRADIENT_MAT<T_MAT,MAT>& v,const T_MAT2& m)
{GRADIENT<typename T_MAT::SCALAR,decltype(DOUBLE_CONTRACT_00_11_VT_M::Type(v.x,m))> r;DOUBLE_CONTRACT_00_11_VT_M()(r.x,v.x,m);return r;}

template<int i,class T_MAT,class A,int d,class VEC> inline void
Extract(VECTOR<A,d>& dx,const GRADIENT_MAT<T_MAT,VEC>& v)
{Extract<i>(dx,v.x);}

template<int i,class T_MAT,class OUT,class VEC>
void Get(OUT& o,const GRADIENT_MAT<T_MAT,VEC>& g)
{
    GET_VEC_HELPER<i>::f(o,g.x);
}

}
}
#endif
