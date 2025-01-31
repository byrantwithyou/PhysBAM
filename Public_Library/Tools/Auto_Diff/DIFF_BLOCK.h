//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIFF_BLOCK
//##################################################################### 
#ifndef __DIFF_BLOCK__
#define __DIFF_BLOCK__
#include <Core/Math_Tools/FIXED_NUMBER.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
#include <Tools/Tensors/DIAGONAL_TENSOR.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <type_traits>
using std::enable_if;

namespace PhysBAM{
template<class T> struct IS_SCALAR{static const bool value=is_scalar<T>::value;};
template<class T,int m> struct IS_SCALAR<FIXED_NUMBER<T,m> >{static const bool value=true;};
#define DECLARE_BLOCK(name,check_type) \
template<class OBJ> \
struct BLOCK_##name \
{ \
    static_assert(IS_##check_type<OBJ>::value,"BLOCK_"#name" must be templatized over a "#check_type); \
    const static bool block_tag=true; \
    OBJ obj; \
    BLOCK_##name() :obj() {}                          \
    explicit BLOCK_##name(const OBJ& obj) :obj(obj){} \
}; \
template<class OBJ> BLOCK_##name<OBJ> Make_##name(const OBJ& obj){return BLOCK_##name<OBJ>(obj);};

DECLARE_BLOCK(sS,SCALAR);
DECLARE_BLOCK(sV,VECTOR);
DECLARE_BLOCK(vS,VECTOR);
DECLARE_BLOCK(vV,MATRIX);
DECLARE_BLOCK(mS,MATRIX);
DECLARE_BLOCK(mV,TENSOR);

DECLARE_BLOCK(sSS,SCALAR);
DECLARE_BLOCK(sVS,VECTOR);
DECLARE_BLOCK(vSS,VECTOR);
DECLARE_BLOCK(vVS,MATRIX);
DECLARE_BLOCK(sSV,VECTOR);
DECLARE_BLOCK(sVV,MATRIX);
DECLARE_BLOCK(vSV,MATRIX);
DECLARE_BLOCK(vVV,TENSOR);

template<template<class OBJ> class BLOCK,class OBJ0,class OBJ1> auto
Add_BB(const BLOCK<OBJ0>& b0,const BLOCK<OBJ1>& b1,enable_if_t<BLOCK<OBJ0>::block_tag,void*> =0)
{
    return BLOCK<typename remove_const<typename remove_reference<decltype(b0.obj+b1.obj)>::type>::type>(b0.obj+b1.obj);
}

template<template<class OBJ> class BLOCK,class OBJ0,class OBJ1> auto
Sub_BB(const BLOCK<OBJ0>& b0,const BLOCK<OBJ1>& b1,enable_if_t<BLOCK<OBJ0>::block_tag,void*> =0)
{
    return BLOCK<typename remove_const<typename remove_reference<decltype(b0.obj-b1.obj)>::type>::type>(b0.obj-b1.obj);
}

template<template<class OBJ> class BLOCK,class OBJ> auto
Neg_B(const BLOCK<OBJ>& b,enable_if_t<BLOCK<OBJ>::block_tag,void*> =0)
{
    return BLOCK<typename remove_const<typename remove_reference<decltype(-b.obj)>::type>::type>(-b.obj);
}

template<template<class OBJ> class BLOCK,class OBJ,class T> auto
Mul_BS(const BLOCK<OBJ>& b,const T& s,enable_if_t<BLOCK<OBJ>::block_tag && IS_SCALAR<T>::value,void*> =0)
{
    return BLOCK<typename remove_const<typename remove_reference<decltype(b.obj*s)>::type>::type>(b.obj*s);
}

template<template<class OBJ> class BLOCK,class OBJ,class T> auto
Div_BS(const BLOCK<OBJ>& b,const T& s,enable_if_t<BLOCK<OBJ>::block_tag && IS_SCALAR<T>::value,void*> =0)
{
    return BLOCK<typename remove_const<typename remove_reference<decltype(b.obj/s)>::type>::type>(b.obj/s);
}

template<class OBJ,class MAT> auto
Mul_MB(const MAT& mat,const BLOCK_vV<OBJ>& block,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vV(mat*block.obj);
}

template<class OBJ,class MAT> auto
Mul_MB(const MAT& mat,const BLOCK_vS<OBJ>& block,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vS(mat*block.obj);
}

template<class OBJ,class VEC> auto
Contract_0(const BLOCK_vSS<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_sSS(block.obj.Dot(vec));
}

template<class OBJ,class VEC> auto
Contract_0(const BLOCK_vVS<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_sVS(Transpose_Times(block.obj,vec));
}

template<class OBJ,class VEC> auto
Contract_0(const BLOCK_vSV<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_sSV(Transpose_Times(block.obj,vec));
}

template<class OBJ,class VEC> auto
Contract_0(const BLOCK_vVV<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_sVV(Contract<0>(block.obj,vec));
}

template<class OBJ,class MAT> auto
Contract_00(const BLOCK_vSS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vSS(Transpose_Times(mat,block.obj));
}

template<class OBJ,class MAT> auto
Contract_00(const BLOCK_vVS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vVS(Transpose_Times(mat,block.obj));
}

template<class OBJ,class MAT> auto
Contract_00(const BLOCK_vSV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vSV(Transpose_Times(mat,block.obj));
}

template<class OBJ,class MAT> auto
Contract_00(const BLOCK_vVV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vVV(Contract<0,0>(block.obj,mat));
}

template<class OBJ> auto
Transpose_Times_Self(const BLOCK_vS<OBJ>& block)
{
    return Make_sSS(block.obj.Dot(block.obj));
}

template<class OBJ> auto
Transpose_Times_Self(const BLOCK_vV<OBJ>& block)
{
    return Make_sVV(Transpose_Times_Self(block.obj));
}

template<class OBJ,class VEC> auto
Transpose_Times(const BLOCK_vS<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_sS(block.obj.Dot(vec));
}

template<class OBJ,class VEC> auto
Transpose_Times(const BLOCK_vV<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_sV(Transpose_Times(block.obj,vec));
}

template<class OBJ,class VEC> auto
Tensor_Product_0(const BLOCK_sVV<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_vVV(Tensor_Product<0>(block.obj,vec));
}

template<class OBJ,class VEC> auto
Tensor_Product_0(const BLOCK_sVS<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_vVS(Outer_Product(vec,block.obj));
}

template<class OBJ,class VEC> auto
Tensor_Product_0(const BLOCK_sSV<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_vSV(Outer_Product(vec,block.obj));
}

template<class OBJ,class VEC> auto
Tensor_Product_0(const BLOCK_sSS<OBJ>& block,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_vSS(vec*block.obj);
}

template<class OBJ0,class VEC> auto
Tensor_Product_1(const BLOCK_vV<OBJ0>& block0,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_mV(Tensor_Product<1>(block0.obj,vec));
}

template<class OBJ0,class VEC> auto
Tensor_Product_1(const BLOCK_vS<OBJ0>& block0,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_mS(Outer_Product(block0.obj,vec));
}

template<class OBJ0,class VEC> auto
Tensor_Product_0(const BLOCK_vV<OBJ0>& block0,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_mV(Tensor_Product<0>(block0.obj,vec));
}

template<class OBJ0,class VEC> auto
Tensor_Product_0(const BLOCK_vS<OBJ0>& block0,const VEC& vec,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_mS(Outer_Product(vec,block0.obj));
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_1(const BLOCK_vV<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_vVV(Tensor_Product<1>(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_1(const BLOCK_vV<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    return Make_vSV(block0.obj*block1.obj);
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_1(const BLOCK_vS<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_vVS(Outer_Product(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_1(const BLOCK_vS<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    return Make_vSS(block0.obj*block1.obj);
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_2(const BLOCK_vV<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_vVV(Tensor_Product<2>(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_2(const BLOCK_vV<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    return Make_vVS(block0.obj*block1.obj);
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_2(const BLOCK_vS<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_vSV(Outer_Product(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Tensor_Product_2(const BLOCK_vS<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    return Make_vSS(block0.obj*block1.obj);
}

template<class OBJ,class MAT> auto
Tensor_Product_2(const MAT& mat,const BLOCK_sV<OBJ>& block,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mV(Tensor_Product<2>(mat,block.obj));
}

template<class OBJ,class MAT> auto
Tensor_Product_2(const MAT& mat,const BLOCK_sS<OBJ>& block,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mS(mat*block.obj);
}

template<class OBJ> auto
Outer_Product(const BLOCK_sV<OBJ>& block)
{
    return Make_sVV(Outer_Product(block.obj));
}

template<class OBJ> auto
Outer_Product(const BLOCK_sS<OBJ>& block)
{
    return Make_sSS(block.obj*block.obj);
}

template<class OBJ,class VEC> auto
Outer_Product(const VEC& vec,const BLOCK_sV<OBJ>& block,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_vV(Outer_Product(vec,block.obj));
}

template<class OBJ,class VEC> auto
Outer_Product(const VEC& vec,const BLOCK_sS<OBJ>& block,enable_if_t<IS_VECTOR<VEC>::value,void*> =0)
{
    return Make_vS(block.obj*vec);
}

template<class OBJ0,class OBJ1> auto
Outer_Product(const BLOCK_sV<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_sVV(Outer_Product(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Outer_Product(const BLOCK_sV<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    return Make_sVS(block0.obj*block1.obj);
}

template<class OBJ0,class OBJ1> auto
Outer_Product(const BLOCK_sS<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_sSV(block0.obj*block1.obj);
}

template<class OBJ0,class OBJ1> auto
Outer_Product(const BLOCK_sS<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    return Make_sSS(block0.obj*block1.obj);
}

template<class OBJ0,class OBJ1> auto
Symmetric_Outer_Product(const BLOCK_sV<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_sVV(Symmetric_Outer_Product(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Symmetric_Outer_Product(const BLOCK_sS<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    auto a=block0.obj*block1.obj;
    return Make_sSS(a+a);
}

template<class OBJ,class MAT> auto
Contract_01(const BLOCK_vVV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vVV(Contract<0,1>(block.obj,mat));
}

template<class OBJ,class MAT> auto
Contract_01(const BLOCK_vVS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vVS(mat*block.obj);
}

template<class OBJ,class MAT> auto
Contract_01(const BLOCK_vSV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vSV(mat*block.obj);
}

template<class OBJ,class MAT> auto
Contract_01(const BLOCK_vSS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_vSS(mat*block.obj);
}

template<class OBJ0,class OBJ1> auto
Transpose_Times(const BLOCK_vV<OBJ0>& block0,const BLOCK_vV<OBJ1>& block1)
{
    return Make_sVV(Transpose_Times(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Transpose_Times(const BLOCK_vS<OBJ0>& block0,const BLOCK_vV<OBJ1>& block1)
{
    return Make_sSV(Transpose_Times(block1.obj,block0.obj));
}

template<class OBJ0,class OBJ1> auto
Transpose_Times(const BLOCK_vV<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1)
{
    return Make_sVS(Transpose_Times(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Transpose_Times(const BLOCK_vS<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1)
{
    return Make_sSS(block0.obj.Dot(block1.obj));
}

template<class OBJ0,class OBJ1> auto
Symmetric_Transpose_Times(const BLOCK_vV<OBJ0>& block0,const BLOCK_vV<OBJ1>& block1)
{
    return Make_sVV(Symmetric_Transpose_Times(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Symmetric_Transpose_Times(const BLOCK_vS<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1)
{
    auto a=block0.obj.Dot(block1.obj);
    return Make_sSS(a+a);
}

template<class OBJ0,class OBJ1> auto
Symmetric_Tensor_Product_12(const BLOCK_vV<OBJ0>& block0,const BLOCK_sV<OBJ1>& block1)
{
    return Make_vVV(Symmetric_Tensor_Product<1,2>(block0.obj,block1.obj));
}

template<class OBJ0,class OBJ1> auto
Symmetric_Tensor_Product_12(const BLOCK_vS<OBJ0>& block0,const BLOCK_sS<OBJ1>& block1)
{
    auto a=block0.obj*block1.obj;
    return Make_vSS(a+a);
}

template<class TEN,class OBJ0,class OBJ1> auto
Double_Contract_12(const TEN& a,const BLOCK_vV<OBJ0>& block0,const BLOCK_vV<OBJ1>& block1,
    enable_if_t<IS_TENSOR<TEN>::value,void*> =0)
{
    return Make_vVV(Contract<2,0>(Contract<1,0>(a,block0.obj),block1.obj));
}

template<class TEN,class OBJ0,class OBJ1> auto
Double_Contract_12(const TEN& a,const BLOCK_vV<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1,
    enable_if_t<IS_TENSOR<TEN>::value,void*> =0)
{
    return Make_vVS(Contract<2>(a,block1.obj)*block0.obj);
}

template<class TEN,class OBJ0,class OBJ1> auto
Double_Contract_12(const TEN& a,const BLOCK_vS<OBJ0>& block0,const BLOCK_vV<OBJ1>& block1,
    enable_if_t<IS_TENSOR<TEN>::value,void*> =0)
{
    return Make_vSV(Contract<1>(a,block0.obj)*block1.obj);
}

template<class T,class OBJ0,class OBJ1> auto
Double_Contract_12(const PERMUTATION_TENSOR<T>& a,const BLOCK_vS<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1)
{
    return Make_vSS(a.x*Cross_Product(block0.obj,block1.obj));
}

template<class T,int d,class OBJ0,class OBJ1> auto
Double_Contract_12(const DIAGONAL_TENSOR<T,d>& a,const BLOCK_vS<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1)
{
    return Make_vSS(a.v*block0.obj*block1.obj);
}

template<class OBJ0,class OBJ1,class T_TEN> auto
Symmetric_Double_Contract_12(const T_TEN& a,const BLOCK_vV<OBJ0>& block0,const BLOCK_vV<OBJ1>& block1)
{
    return Make_vVV(Symmetric_Double_Contract_12(a,block0.obj,block1.obj));
}

template<class T,class OBJ0,class OBJ1> auto
Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& a,const BLOCK_vS<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1)
{
    return Make_vSS(a.x*2*Cross_Product(block0.obj,block1.obj));
}

template<class T,int m,class OBJ0,class OBJ1> auto
Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& a,const BLOCK_vS<OBJ0>& block0,const BLOCK_vS<OBJ1>& block1)
{
    return Make_vSS(a.v*2*block0.obj*block1.obj);
}

template<class OBJ0> auto
Twice_Symmetric_Part_01(const BLOCK_mV<OBJ0>& block0)
{
    return Make_mV(Twice_Symmetric_Part<0,1>(block0.obj));
}

template<class OBJ0> auto
Twice_Symmetric_Part_01(const BLOCK_mS<OBJ0>& block0)
{
    return Make_mS(block0.obj.Twice_Symmetric_Part());
}

template<class OBJ0> auto
Contract_01(const BLOCK_mV<OBJ0>& block0)
{
    return Make_sV(Contract<0,1>(block0.obj));
}

template<class OBJ0> auto
Contract_01(const BLOCK_mS<OBJ0>& block0)
{
    return Make_sS(block0.obj.Trace());
}

template<class OBJ0> auto
Transposed_01(const BLOCK_mV<OBJ0>& block0)
{
    return Make_mV(Transposed<0,1>(block0.obj));
}

template<class OBJ0> auto
Transposed_01(const BLOCK_mS<OBJ0>& block0)
{
    return Make_mS(block0.obj.Transposed());
}

template<class OBJ,class MAT> auto
Contract_00(const BLOCK_mV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mV(Contract<0,0>(block.obj,mat));
}

template<class OBJ,class MAT> auto
Contract_00(const BLOCK_mS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mS(Transpose_Times(block.obj,mat));
}

template<class OBJ,class MAT> auto
Contract_01(const BLOCK_mV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mV(Contract<0,1>(block.obj,mat));
}

template<class OBJ,class MAT> auto
Contract_01(const BLOCK_mS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mS((mat*block.obj).Transposed());
}

template<class OBJ,class MAT> auto
Contract_10(const BLOCK_mV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mV(Contract<1,0>(block.obj,mat));
}

template<class OBJ,class MAT> auto
Contract_10(const BLOCK_mS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mS(block.obj*mat);
}

template<class OBJ,class MAT> auto
Contract_11(const BLOCK_mV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mV(Contract<1,1>(block.obj,mat));
}

template<class OBJ,class MAT> auto
Contract_11(const BLOCK_mS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_mS(Times_Transpose(block.obj,mat));
}

template<class OBJ,class MAT> auto
Double_Contract_00_11(const BLOCK_mV<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_sV(Contract<0,1>(Contract<1,1>(block.obj,mat)));
}

template<class OBJ,class MAT> auto
Double_Contract_00_11(const BLOCK_mS<OBJ>& block,const MAT& mat,enable_if_t<IS_MATRIX<MAT>::value,void*> =0)
{
    return Make_sS(block.obj.Double_Contract(mat));
}

template<template<class OBJ> class BLOCK,class OBJ0,class OBJ1> auto
Choose(const BLOCK<OBJ0>& block0,const BLOCK<OBJ1>& block1)
    -> BLOCK<decltype(Choose(block0.obj,block1.obj))>;

template<template<class OBJ> class BLOCK,class OBJ> auto
Choose_Zero(const BLOCK<OBJ>& block)
    -> BLOCK<decltype(Choose_Zero(block.obj))>;

template<int k,class T,int... dims> struct ZERO_BLOCK_TYPE;
template<class T> struct ZERO_BLOCK_TYPE<1,T,-1,-1> {typedef BLOCK_sS<FIXED_NUMBER<T,0> > TYPE;};
template<class T,int n> struct ZERO_BLOCK_TYPE<1,T,-1,n> {typedef BLOCK_sV<ZERO_VECTOR<T,n> > TYPE;};
template<class T,int m> struct ZERO_BLOCK_TYPE<1,T,m,-1> {typedef BLOCK_vS<ZERO_VECTOR<T,m> > TYPE;};
template<class T,int m,int n> struct ZERO_BLOCK_TYPE<1,T,m,n> {typedef BLOCK_vV<ZERO_MATRIX<T,m,n> > TYPE;};
template<class T,int m,int n> struct ZERO_BLOCK_TYPE<1,T,m,n,-1> {typedef BLOCK_mS<ZERO_MATRIX<T,m,n> > TYPE;};
template<class T,int m,int n,int p> struct ZERO_BLOCK_TYPE<1,T,m,n,p> {typedef BLOCK_mV<ZERO_TENSOR<T,m,n,p> > TYPE;};

template<class T> struct ZERO_BLOCK_TYPE<2,T,-1,-1,-1> {typedef BLOCK_sSS<FIXED_NUMBER<T,0> > TYPE;};
template<class T,int n> struct ZERO_BLOCK_TYPE<2,T,-1,n,-1> {typedef BLOCK_sVS<ZERO_VECTOR<T,n> > TYPE;};
template<class T,int m> struct ZERO_BLOCK_TYPE<2,T,m,-1,-1> {typedef BLOCK_vSS<ZERO_VECTOR<T,m> > TYPE;};
template<class T,int m,int n> struct ZERO_BLOCK_TYPE<2,T,m,n,-1> {typedef BLOCK_vVS<ZERO_MATRIX<T,m,n> > TYPE;};

template<class T,int p> struct ZERO_BLOCK_TYPE<2,T,-1,-1,p> {typedef BLOCK_sSV<ZERO_VECTOR<T,p> > TYPE;};
template<class T,int n,int p> struct ZERO_BLOCK_TYPE<2,T,-1,n,p> {typedef BLOCK_sVV<ZERO_MATRIX<T,n,p> > TYPE;};
template<class T,int m,int p> struct ZERO_BLOCK_TYPE<2,T,m,-1,p> {typedef BLOCK_vSV<ZERO_MATRIX<T,m,p> > TYPE;};
template<class T,int m,int n,int p> struct ZERO_BLOCK_TYPE<2,T,m,n,p> {typedef BLOCK_vVV<ZERO_TENSOR<T,m,n,p> > TYPE;};

template<int k,class T,int... dims> struct IDENTITY_BLOCK_TYPE;
template<class T> struct IDENTITY_BLOCK_TYPE<1,T,-1,-1> {typedef BLOCK_sS<FIXED_NUMBER<T,1> > TYPE;};
template<class T,int m> struct IDENTITY_BLOCK_TYPE<1,T,m,m> {typedef BLOCK_vV<IDENTITY_MATRIX<T,m> > TYPE;};

template<class T> struct IDENTITY_BLOCK_TYPE<2,T,-1,-1,-1> {typedef BLOCK_sSS<FIXED_NUMBER<T,1> > TYPE;};
template<class T,int n> struct IDENTITY_BLOCK_TYPE<2,T,-1,n,n> {typedef BLOCK_sVV<IDENTITY_MATRIX<T,n> > TYPE;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES{static const bool value=false;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sS<A>,BLOCK_sS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sV<A>,BLOCK_sS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sS<A>,BLOCK_sV<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sV<A>,BLOCK_sV<B> >{static const bool value=true;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vS<A>,BLOCK_vS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vV<A>,BLOCK_vS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vS<A>,BLOCK_vV<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vV<A>,BLOCK_vV<B> >{static const bool value=true;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_mS<A>,BLOCK_mS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_mV<A>,BLOCK_mS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_mS<A>,BLOCK_mV<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_mV<A>,BLOCK_mV<B> >{static const bool value=true;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sSS<A>,BLOCK_sSS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sSV<A>,BLOCK_sSS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sSS<A>,BLOCK_sSV<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sSV<A>,BLOCK_sSV<B> >{static const bool value=true;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sVS<A>,BLOCK_sVS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sVV<A>,BLOCK_sVS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sVS<A>,BLOCK_sVV<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_sVV<A>,BLOCK_sVV<B> >{static const bool value=true;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vSS<A>,BLOCK_vSS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vSV<A>,BLOCK_vSS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vSS<A>,BLOCK_vSV<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vSV<A>,BLOCK_vSV<B> >{static const bool value=true;};

template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vVS<A>,BLOCK_vVS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vVV<A>,BLOCK_vVS<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vVS<A>,BLOCK_vVV<B> >{static const bool value=true;};
template<class A,class B> struct ASSERT_SAME_BLOCK_TYPES<BLOCK_vVV<A>,BLOCK_vVV<B> >{static const bool value=true;};

}
#endif
