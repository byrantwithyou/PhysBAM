//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VEC_HOLDER
//##################################################################### 
#ifndef __VEC_HOLDER__
#define __VEC_HOLDER__

#include <Tools/Auto_Diff/DIFF_BLOCK.h>
#include <Tools/Auto_Diff/DIFF_LAYOUT.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class OP,class ...Args> struct RET {typedef decltype(OP()(typename remove_reference<Args>::type()...)) TYPE;};

template<int n> struct ARG {template<class T,class ...Args> auto operator()(const T& t,Args&&... args)->decltype(ARG<n-1>()(args...)){return ARG<n-1>()(args...);}};
template<> struct ARG<0> {template<class T,class ...Args> const T&operator()(const T& t,Args&&... args){return t;}};

#define MK_FUN_1(NAME,fun) template<class A> struct NAME {template<class ...Args> auto operator()(Args&&... args) const -> decltype(fun(A()(args...))) {return fun(A()(args...));}}
#define MK_OP_2(NAME,op) template<class A,class B> struct NAME {template<class ...Args> auto operator()(Args&&... args) const -> decltype(A()(args...) op B()(args...)) {return A()(args...) op B()(args...);}}
#define MK_FUN_2(NAME,fun) template<class A,class B> struct NAME {template<class ...Args> auto operator()(Args&&... args) const -> decltype(fun(A()(args...),B()(args...))) {return fun(A()(args...),B()(args...));}}
#define MK_FUN_3(NAME,fun) template<class A,class B,class C> struct NAME {template<class ...Args> auto operator()(Args&&... args) const -> decltype(fun(A()(args...),B()(args...),C()(args...))) {return fun(A()(args...),B()(args...),C()(args...));}}

MK_FUN_1(NEG,-);
MK_OP_2(ADD,+);
MK_OP_2(SUB,-);
MK_OP_2(MUL,*);
MK_OP_2(DIV,/);
MK_FUN_2(TRANS_MUL,Transpose_Times);

struct VEC_END {};

template<class OBJ_IN,class BASE_IN>
struct VEC_HOLDER
{
//    typedef typename OBJ_IN::SCALAR T;
    typedef OBJ_IN OBJ;
    typedef BASE_IN BASE;
    OBJ x;
    BASE z;
};

inline void Fill_From(VEC_END& o,const VEC_END& v) {}
template<class OBJ,class BASE,class OBJ1,class BASE1>
void Fill_From(VEC_HOLDER<OBJ1,BASE1>& o,const VEC_HOLDER<OBJ,BASE>& v)
{
    Fill_From(o.x.obj,v.x.obj);
    Fill_From(o.z,v.z);
}

template<class OUT> void Get(OUT& o,const VEC_END& v,int i) {PHYSBAM_FATAL_ERROR();}
template<class OUT,class OBJ,class BASE>
void Get(OUT& o,const VEC_HOLDER<OBJ,BASE>& v,int i)
{
    if(i>0) return Get(o,v.z,i-1);
    Fill_From(o,v.x.obj);
}

template<int d,int i> struct EXTRACT_VEC_HELPER
{
    template<class TV,class OBJ,class BASE>
    static void f(VECTOR<TV,d>& dx,const VEC_HOLDER<OBJ,BASE>& v)
    {Fill_From(dx(i),v.x.obj);EXTRACT_VEC_HELPER<d,i+1>::f(dx,v.z);}
};

template<int d> struct EXTRACT_VEC_HELPER<d,d>
{
    template<class TV>
    static void f(VECTOR<TV,d>& dx,const VEC_END& v)
    {}
};

template<class TV,int d,class OBJ,class BASE> inline void
Extract(VECTOR<TV,d>& dx,const VEC_HOLDER<OBJ,BASE>& v)
{EXTRACT_VEC_HELPER<d,0>::f(dx,v);}

template<int n,class IN,class OBJ,class BASE> void Set(VEC_HOLDER<OBJ,BASE>& out,const IN& in) {Set_Helper(out,in,(VECTOR<int,n>*)0);}
template<int n,class IN,class OBJ,class BASE> void Set_Helper(VEC_HOLDER<OBJ,BASE>& out,const IN& in,VECTOR<int,n>*) {Set_Helper(out.z,in,(VECTOR<int,n-1>*)0);}
template<class IN,class OBJ,class BASE> void Set_Helper(VEC_HOLDER<OBJ,BASE>& out,const IN& in,VECTOR<int,0>*){out.x=in;}

template<class LAYOUT,int... lead_dims> struct EMPTY_VEC;
template<class T,int d,int... dims,int... lead_dims> struct EMPTY_VEC<DIFF_LAYOUT<T,d,dims...>,lead_dims...>
{
    typedef VEC_HOLDER<typename ZERO_BLOCK_TYPE<T,lead_dims...,d>::TYPE,typename EMPTY_VEC<DIFF_LAYOUT<T,dims...>,lead_dims...>::TYPE> TYPE;
};
template<class T,int... lead_dims> struct EMPTY_VEC<DIFF_LAYOUT<T>,lead_dims...>
{
    typedef VEC_END TYPE;
};

template<int i,class LAYOUT,int... lead_dims> struct ONE_NONZERO_VECTOR;
template<int i,class T,int d,int... dims,int... lead_dims> struct ONE_NONZERO_VECTOR<i,DIFF_LAYOUT<T,d,dims...>,lead_dims...>
{
    typedef VEC_HOLDER<typename ZERO_BLOCK_TYPE<T,lead_dims...,d>::TYPE,typename ONE_NONZERO_VECTOR<i-1,DIFF_LAYOUT<T,dims...>,lead_dims...>::TYPE> TYPE;
};
template<class T,int d,int... dims,int... lead_dims> struct ONE_NONZERO_VECTOR<0,DIFF_LAYOUT<T,d,dims...>,lead_dims...>
{
    typedef VEC_HOLDER<typename IDENTITY_BLOCK_TYPE<T,lead_dims...,d>::TYPE,typename EMPTY_VEC<DIFF_LAYOUT<T,dims...>,lead_dims...>::TYPE> TYPE;
};

template<class OP,class VEC> struct TYPE_VEC_MAP_1;

template<class OP>
struct TYPE_VEC_MAP_1<OP,VEC_END>
{
    template<class ...Args> static VEC_END Type(Args...);
    template<class ...Args> static void Type_Debug(Args...){}
};

template<class OP,class OBJ,class BASE>
struct TYPE_VEC_MAP_1<OP,VEC_HOLDER<OBJ,BASE> >
{
    template<class ...Args> static VEC_HOLDER<decltype(OP()(OBJ(),typename remove_reference<Args>::type()...)),decltype(TYPE_VEC_MAP_1<OP,BASE>::Type(typename remove_reference<Args>::type()...))>Type(Args...);
    template<class ...Args> static void Type_Debug(Args...){OP()(OBJ(),typename remove_reference<Args>::type()...);TYPE_VEC_MAP_1<OP,BASE>::Type_Debug(typename remove_reference<Args>::type()...);}
};

template<class OP> struct VEC_MAP_1
{
    template<class OBJ,class BASE,class ...Args> static decltype(TYPE_VEC_MAP_1<OP,VEC_HOLDER<OBJ,BASE> >::Type(typename remove_reference<Args>::type()...))
    Type(VEC_HOLDER<OBJ,BASE>,Args...);
    template<class ...Args> static VEC_END Type(VEC_END,Args...);
    template<class OBJ,class BASE,class ...Args> static void Type_Debug(VEC_HOLDER<OBJ,BASE>,Args...)
    {TYPE_VEC_MAP_1<OP,VEC_HOLDER<OBJ,BASE> >::Type_Debug(typename remove_reference<Args>::type()...);}
    template<class ...Args> static void Type_Debug(VEC_END,Args...){}

    template<class ...Args> void operator()(VEC_END& out,const VEC_END& in,Args&&... args) const {}
    template<class OBJ,class BASE,class ...Args> void operator()(decltype(Type(VEC_HOLDER<OBJ,BASE>(),typename remove_reference<Args>::type()...))& out,const VEC_HOLDER<OBJ,BASE>& in,Args&&... args) const
    {
        out.x=OP()(in.x,args...);
        (*this)(out.z,in.z,args...);
    };
};

template<class OP,class VEC0,class VEC1> struct TYPE_VEC_MAP_2;

template<class OP>
struct TYPE_VEC_MAP_2<OP,VEC_END,VEC_END>
{
    template<class ...Args> static VEC_END Type(Args...);
    template<class ...Args> static void Type_Debug(Args...){}
};

template<class OP,class OBJ0,class BASE0,class OBJ1,class BASE1>
struct TYPE_VEC_MAP_2<OP,VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1> >
{
    template<class ...Args> static VEC_HOLDER<decltype(OP()(OBJ0(),OBJ1(),typename remove_reference<Args>::type()...)),decltype(TYPE_VEC_MAP_2<OP,BASE0,BASE1>::Type(typename remove_reference<Args>::type()...))>Type(Args...);
    template<class ...Args> static void Type_Debug(Args...){OP()(OBJ0(),OBJ1(),typename remove_reference<Args>::type()...);TYPE_VEC_MAP_2<OP,BASE0,BASE1>::Type_Debug(typename remove_reference<Args>::type()...);}
};

template<class OP> struct VEC_MAP_2
{
    template<class OBJ0,class BASE0,class OBJ1,class BASE1,class ...Args> static decltype(TYPE_VEC_MAP_2<OP,VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1> >::Type(typename remove_reference<Args>::type()...))
    Type(VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1>,Args...);
    template<class ...Args> static VEC_END Type(VEC_END,VEC_END,Args...);

    template<class OBJ0,class BASE0,class OBJ1,class BASE1,class ...Args> static void
    Type_Debug(VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1>,Args...)
    {TYPE_VEC_MAP_2<OP,VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1> >::Type_Debug(typename remove_reference<Args>::type()...);}
    template<class ...Args> static void Type_Debug(VEC_END,VEC_END,Args...){}

    template<class ...Args> void operator()(VEC_END& out,const VEC_END& in0,const VEC_END& in1,Args&&... args) const {}
    template<class OBJ0,class BASE0,class OBJ1,class BASE1,class ...Args> void
    operator()(decltype(Type(VEC_HOLDER<OBJ0,BASE0>(),VEC_HOLDER<OBJ1,BASE1>(),typename remove_reference<Args>::type()...))& out,const VEC_HOLDER<OBJ0,BASE0>& in0,const VEC_HOLDER<OBJ1,BASE1>& in1,Args&&... args) const
    {
        out.x=OP()(in0.x,in1.x,args...);
        (*this)(out.z,in0.z,in1.z,args...);
    };
};

typedef VEC_MAP_1<NEG<ARG<0> > > VEC_NEG;
typedef VEC_MAP_2<ADD<ARG<0>,ARG<1> > > VEC_ADD;
typedef VEC_MAP_2<SUB<ARG<0>,ARG<1> > > VEC_SUB;
typedef VEC_MAP_1<MUL<ARG<0>,ARG<1> > > VEC_SCALE;
typedef VEC_MAP_1<DIV<ARG<0>,ARG<1> > > VEC_SCALE_DIV;
typedef VEC_MAP_1<TRANS_MUL<ARG<0>,ARG<1> > > VEC_TRANSPOSE_TIMES;

MK_FUN_2(TENSOR_PRODUCT_0,Tensor_Product_0);
MK_FUN_2(TENSOR_PRODUCT_1,Tensor_Product_1);
MK_FUN_2(TENSOR_PRODUCT_2,Tensor_Product_2);
MK_FUN_2(CONTRACT_0,Contract_0);
MK_FUN_2(CONTRACT_1,Contract_1);
MK_FUN_2(CONTRACT_2,Contract_2);
MK_FUN_2(CONTRACT_00,Contract_00);
MK_FUN_2(CONTRACT_10,Contract_10);
MK_FUN_2(CONTRACT_01,Contract_01);
MK_FUN_2(CONTRACT_20,Contract_20);
MK_FUN_2(OUTER_PRODUCT,Outer_Product);
MK_FUN_2(CHOOSE,Choose);

typedef VEC_MAP_1<TENSOR_PRODUCT_0<ARG<1>,ARG<0> > > TENSOR_PRODUCT_0_VV_M;
typedef VEC_MAP_1<TENSOR_PRODUCT_0<ARG<0>,ARG<1> > > TENSOR_PRODUCT_0_VM_V;

typedef VEC_MAP_1<TENSOR_PRODUCT_1<ARG<1>,ARG<0> > > TENSOR_PRODUCT_1_VV_M;
typedef VEC_MAP_1<TENSOR_PRODUCT_1<ARG<0>,ARG<1> > > TENSOR_PRODUCT_1_VM_V;

typedef VEC_MAP_1<TENSOR_PRODUCT_2<ARG<1>,ARG<0> > > TENSOR_PRODUCT_2_VV_M;
typedef VEC_MAP_1<TENSOR_PRODUCT_2<ARG<0>,ARG<1> > > TENSOR_PRODUCT_2_VM_V;

typedef VEC_MAP_1<CONTRACT_0<ARG<0>,ARG<1> > > CONTRACT_0_VT_V;
typedef VEC_MAP_1<CONTRACT_0<ARG<1>,ARG<0> > > CONTRACT_0_VV_T;
typedef VEC_MAP_1<CONTRACT_00<ARG<0>,ARG<1> > > CONTRACT_0_VT_M;
typedef VEC_MAP_1<CONTRACT_00<ARG<1>,ARG<0> > > CONTRACT_0_VM_T;
typedef VEC_MAP_1<CONTRACT_01<ARG<1>,ARG<0> > > TRANSPOSE_CONTRACT_0_VM_T;
typedef VEC_MAP_1<CONTRACT_01<ARG<0>,ARG<1> > > VEC_SCALE_REV;

typedef VEC_MAP_1<OUTER_PRODUCT<ARG<0>,ARG<1> > > VEC_OUTER_PRODUCT;
typedef VEC_MAP_1<OUTER_PRODUCT<ARG<1>,ARG<0> > > VEC_OUTER_PRODUCT_REV;
typedef VEC_MAP_2<CHOOSE<ARG<0>,ARG<1> > > VEC_CHOOSE;
}
}
#endif
