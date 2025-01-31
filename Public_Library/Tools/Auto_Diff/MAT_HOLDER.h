//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MAT_HOLDER
//##################################################################### 
#ifndef __MAT_HOLDER__
#define __MAT_HOLDER__

#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/VEC_HOLDER.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

struct MAT_END {};

template<class OBJ_IN,class COL_IN,class BASE_IN>
struct MAT_HOLDER
{
//    typedef typename OBJ_IN::SCALAR T;
    typedef OBJ_IN OBJ;
    typedef COL_IN COL;
    typedef BASE_IN BASE;
    OBJ x;
    COL y;
    BASE z;
};

template<class V>
struct ASSERT_VALID_BLOCK_TYPES_MAT
{
    static const bool value=ASSERT_SAME_BLOCK_TYPES_VEC<typename V::OBJ,typename V::COL>::value&&ASSERT_VALID_BLOCK_TYPES_MAT<typename V::BASE>::value;
};

template<>
struct ASSERT_VALID_BLOCK_TYPES_MAT<MAT_END>
{
    static const bool value=true;
};

inline void Fill_From(MAT_END& o,const MAT_END& m) {}
template<class OBJ,class COL,class BASE,class OBJ1,class COL1,class BASE1>
void Fill_From(MAT_HOLDER<OBJ1,COL1,BASE1>& o,const MAT_HOLDER<OBJ,COL,BASE>& m)
{
    Fill_From(o.x.obj,m.x.obj);
    Fill_From(o.y,m.y);
    Fill_From(o.z,m.z);
}

template<int i,int j>
struct GET_MAT_HELPER
{
    template<class OUT,class OBJ,class COL,class BASE>
    static void f(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& h){GET_MAT_HELPER<i-1,j-1>::f(o,h.z);}
};

template<int i>
struct GET_MAT_HELPER<i,0>
{
    template<class OUT,class OBJ,class COL,class BASE>
    static void f(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& h){GET_VEC_HELPER<i-1>::ft(o,h.y);}
};

template<int j>
struct GET_MAT_HELPER<0,j>
{
    template<class OUT,class OBJ,class COL,class BASE>
    static void f(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& h){GET_VEC_HELPER<j-1>::f(o,h.y);}
};

template<>
struct GET_MAT_HELPER<0,0>
{
    template<class OUT,class OBJ,class COL,class BASE>
    static void f(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& h){Fill_From(o,h.x.obj);}
};

template<int i,int j,class OUT,class OBJ,class COL,class BASE>
void Get(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& v)
{
    GET_MAT_HELPER<i,j>(o,v);
}

template<int i,int j,int ni,int nj,class A,class OBJ>
inline enable_if_t<(i>=0 && i<ni && j>=0 && j<nj)>
Extract_Entry_Helper(MATRIX<A,ni,nj>& ddx,const OBJ& h)
{
    Fill_From(ddx(i,j),h.obj);
}
template<int i,int j,int ni,int nj,class A,class OBJ>
inline enable_if_t<!(i>=0 && i<ni && j>=0 && j<nj)>
Extract_Entry_Helper(MATRIX<A,ni,nj>& ddx,const OBJ& h) {}

template<int i,int j,int ni,int nj,class A,class OBJ>
inline enable_if_t<(i>=0 && i<ni && j>=0 && j<nj)>
Extract_Entry_Transpose_Helper(MATRIX<A,ni,nj>& ddx,const OBJ& h)
{
    Fill_From_Transpose(ddx(i,j),h.obj);
}

template<int i,int j,int ni,int nj,class A,class OBJ>
inline enable_if_t<!(i>=0 && i<ni && j>=0 && j<nj)>
Extract_Entry_Transpose_Helper(MATRIX<A,ni,nj>& ddx,const OBJ& h) {}

template<int i,int j,int ni,int nj,class A,class OBJ,class BASE> inline void
Extract_Matrix_Row_Helper(MATRIX<A,ni,nj>& ddx,const VEC_HOLDER<OBJ,BASE>& h)
{
    Extract_Entry_Helper<i,j>(ddx,h.x);
    Extract_Matrix_Row_Helper<i,j+1>(ddx,h.z);
}

template<int i,int j,int ni,int nj,class A> inline void
Extract_Matrix_Row_Helper(MATRIX<A,ni,nj>& ddx,const VEC_END& h) {}

template<int i,int j,int ni,int nj,class A,class OBJ,class BASE> inline void
Extract_Matrix_Column_Helper(MATRIX<A,ni,nj>& ddx,const VEC_HOLDER<OBJ,BASE>& h)
{
    Extract_Entry_Transpose_Helper<i,j>(ddx,h.x);
    Extract_Matrix_Column_Helper<i+1,j>(ddx,h.z);
}

template<int i,int j,int ni,int nj,class A> inline void
Extract_Matrix_Column_Helper(MATRIX<A,ni,nj>& ddx,const VEC_END& h) {}

template<int i,int j,int ni,int nj,class A,class OBJ,class COL,class BASE> inline void
Extract_Matrix_Helper(MATRIX<A,ni,nj>& ddx,const MAT_HOLDER<OBJ,COL,BASE>& h)
{
    Extract_Entry_Helper<i,j>(ddx,h.x);
    Extract_Matrix_Row_Helper<i,j+1>(ddx,h.y);
    Extract_Matrix_Column_Helper<i+1,j>(ddx,h.y);
    Extract_Matrix_Helper<i+1,j+1>(ddx,h.z);
}

template<int i,int j,int ni,int nj,class A> inline void
Extract_Matrix_Helper(MATRIX<A,ni,nj>& ddx,const MAT_END& h) {}

template<int i,int j,int ni,int nj,class A,class OBJ,class COL,class BASE> inline void
Extract(MATRIX<A,ni,nj>& ddx,const MAT_HOLDER<OBJ,COL,BASE>& h)
{Extract_Matrix_Helper<-i,-j>(ddx,h);}

template<class LAYOUT,int... lead_dims> struct EMPTY_MAT;
template<class T,int d,int... dims,int... lead_dims> struct EMPTY_MAT<DIFF_LAYOUT<T,d,dims...>,lead_dims...>
{
    typedef MAT_HOLDER<
        typename ZERO_BLOCK_TYPE<2,T,lead_dims...,d,d>::TYPE,
        typename EMPTY_VEC<2,DIFF_LAYOUT<T,dims...>,lead_dims...,d>::TYPE,
        typename EMPTY_MAT<DIFF_LAYOUT<T,dims...>,lead_dims...>::TYPE
        > TYPE;
};
template<class T,int... lead_dims> struct EMPTY_MAT<DIFF_LAYOUT<T>,lead_dims...>
{
    typedef MAT_END TYPE;
};

template<class OP,class VEC_OP,class MAT> struct TYPE_MAT_MAP_1;

template<class OP,class VEC_OP>
struct TYPE_MAT_MAP_1<OP,VEC_OP,MAT_END>
{template<class ...Args> static MAT_END Type(Args&&...);};

template<class OP,class VEC_OP,class OBJ,class COL,class BASE>
struct TYPE_MAT_MAP_1<OP,VEC_OP,MAT_HOLDER<OBJ,COL,BASE> >
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ(),typename remove_reference<Args>::type()...)),decltype(VEC_OP::Type(COL(),typename remove_reference<Args>::type()...)),decltype(TYPE_MAT_MAP_1<OP,VEC_OP,BASE>::Type(typename remove_reference<Args>::type()...))> Type(Args&&...);};

template<class OP,class VEC_OP> struct MAT_MAP_1
{
    template<class OBJ,class COL,class BASE,class ...Args>
    static decltype(TYPE_MAT_MAP_1<OP,VEC_OP,MAT_HOLDER<OBJ,COL,BASE> >::Type(typename remove_reference<Args>::type()...))
    Type(const MAT_HOLDER<OBJ,COL,BASE>&,Args&&...);
    template<class ...Args> static MAT_END Type(MAT_END,Args&&...);
    template<class ...Args> static DIFF_UNUSED Type(DIFF_UNUSED,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const MAT_END& in,Args&&... args) const {}
    template<class OBJ,class COL,class BASE,class ...Args> void operator()(decltype(Type(MAT_HOLDER<OBJ,COL,BASE>(),typename remove_reference<Args>::type()...))& out,const MAT_HOLDER<OBJ,COL,BASE>& in,Args&&... args) const
    {
        out.x=OP()(in.x,args...);
        VEC_OP()(out.y,in.y,args...);
        (*this)(out.z,in.z,args...);
    };
};

template<class OP,class VEC_OP,class MAT0,class MAT1> struct TYPE_MAT_MAP_2;

template<class OP,class VEC_OP>
struct TYPE_MAT_MAP_2<OP,VEC_OP,MAT_END,MAT_END>
{template<class ...Args> static MAT_END Type(Args&&...);};

template<class OP,class VEC_OP,class OBJ0,class COL0,class BASE0,class OBJ1,class COL1,class BASE1>
struct TYPE_MAT_MAP_2<OP,VEC_OP,MAT_HOLDER<OBJ0,COL0,BASE0>,MAT_HOLDER<OBJ1,COL1,BASE1> >
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ0(),OBJ1(),typename remove_reference<Args>::type()...)),decltype(VEC_OP::Type(COL0(),COL1(),typename remove_reference<Args>::type()...)),decltype(TYPE_MAT_MAP_2<OP,VEC_OP,BASE0,BASE1>::Type(typename remove_reference<Args>::type()...))>Type(Args&&...);};

template<class OP,class VEC_OP> struct MAT_MAP_2
{
    template<class OBJ0,class COL0,class BASE0,class OBJ1,class COL1,class BASE1,class ...Args>
    static decltype(TYPE_MAT_MAP_2<OP,VEC_OP,MAT_HOLDER<OBJ0,COL0,BASE0>,MAT_HOLDER<OBJ1,COL1,BASE1> >::Type(typename remove_reference<Args>::type()...))
    Type(const MAT_HOLDER<OBJ0,COL0,BASE0>&,const MAT_HOLDER<OBJ1,COL1,BASE1>&,Args&&...);
    template<class ...Args> static MAT_END Type(MAT_END,MAT_END,Args&&...);
    template<class ...Args> static DIFF_UNUSED Type(DIFF_UNUSED,DIFF_UNUSED,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const MAT_END& in0,const MAT_END& in1,Args&&... args) const {}
    template<class OBJ0,class COL0,class BASE0,class OBJ1,class COL1,class BASE1,class ...Args> void
    operator()(decltype(Type(MAT_HOLDER<OBJ0,COL0,BASE0>(),MAT_HOLDER<OBJ1,COL1,BASE1>(),typename remove_reference<Args>::type()...))& out,const MAT_HOLDER<OBJ0,COL0,BASE0>& in0,const MAT_HOLDER<OBJ1,COL1,BASE1>& in1,Args&&... args) const
    {
        out.x=OP()(in0.x,in1.x,args...);
        VEC_OP()(out.y,in0.y,in1.y,args...);
        (*this)(out.z,in0.z,in1.z,args...);
    };
};

template<class OP,class VEC_OP,class VEC0,class VEC1> struct TYPE_SYM_OUTER_MAP;

template<class OP,class VEC_OP>
struct TYPE_SYM_OUTER_MAP<OP,VEC_OP,VEC_END,VEC_END>
{template<class ...Args> static MAT_END Type(Args&&...);};

template<class OP,class VEC_OP,class OBJ0,class BASE0,class OBJ1,class BASE1>
struct TYPE_SYM_OUTER_MAP<OP,VEC_OP,VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1> >
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ0(),OBJ1(),typename remove_reference<Args>::type()...)),decltype(VEC_OP::Type(BASE0(),BASE1(),OBJ0(),OBJ1(),typename remove_reference<Args>::type()...)),decltype(TYPE_SYM_OUTER_MAP<OP,VEC_OP,BASE0,BASE1>::Type(typename remove_reference<Args>::type()...))> Type(Args&&...);};

template<class OP,class VEC_OP> struct SYM_OUTER_MAP
{
    template<class ...Args> static MAT_END Type(VEC_END,VEC_END,Args&&... args);
    template<class OBJ0,class BASE0,class OBJ1,class BASE1,class ...Args> static
    decltype(TYPE_SYM_OUTER_MAP<OP,VEC_OP,VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1> >::Type(typename remove_reference<Args>::type()...)) Type(const VEC_HOLDER<OBJ0,BASE0>&,const VEC_HOLDER<OBJ1,BASE1>&,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const VEC_END& in0,const VEC_END& in1,Args&&... args) const {}
    template<class OBJ0,class BASE0,class OBJ1,class BASE1,class ...Args> void
    operator()(decltype(Type(VEC_HOLDER<OBJ0,BASE0>(),VEC_HOLDER<OBJ1,BASE1>(),typename remove_reference<Args>::type()...))& out,const VEC_HOLDER<OBJ0,BASE0>& in0,const VEC_HOLDER<OBJ1,BASE1>& in1,Args&&... args) const
    {
        out.x=OP()(in0.x,in1.x,args...);
        VEC_OP()(out.y,in0.z,in1.z,in0.x,in1.x,args...);
        (*this)(out.z,in0.z,in1.z,args...);
    };
};

template<class OP,class VEC_OP,class VEC> struct TYPE_OUTER_MAP;

template<class OP,class VEC_OP>
struct TYPE_OUTER_MAP<OP,VEC_OP,VEC_END>
{template<class ...Args> static MAT_END Type(Args&&...);};

template<class OP,class VEC_OP,class OBJ,class BASE>
struct TYPE_OUTER_MAP<OP,VEC_OP,VEC_HOLDER<OBJ,BASE> >
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ(),typename remove_reference<Args>::type()...)),decltype(VEC_OP::Type(BASE(),OBJ(),typename remove_reference<Args>::type()...)),decltype(TYPE_OUTER_MAP<OP,VEC_OP,BASE>::Type(typename remove_reference<Args>::type()...))> Type(Args&&...);};

template<class OP,class VEC_OP> struct OUTER_MAP
{
    template<class ...Args> static MAT_END Type(VEC_END,VEC_END,Args&&... args);
    template<class OBJ,class BASE,class ...Args> static
    decltype(TYPE_OUTER_MAP<OP,VEC_OP,VEC_HOLDER<OBJ,BASE> >::Type(typename remove_reference<Args>::type()...)) Type(const VEC_HOLDER<OBJ,BASE>&,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const VEC_END& in,Args&&... args) const {}
    template<class OBJ,class BASE,class ...Args> void
    operator()(decltype(Type(VEC_HOLDER<OBJ,BASE>(),typename remove_reference<Args>::type()...))& out,const VEC_HOLDER<OBJ,BASE>& in,Args&&... args) const
    {
        out.x=OP()(in.x,args...);
        VEC_OP()(out.y,in.z,in.x,args...);
        (*this)(out.z,in.z,args...);
    };
};

MK_FUN_1(OUTER_PRODUCT_EQ,Outer_Product);
MK_FUN_2(SYMMETRIC_OUTER_PRODUCT,Symmetric_Outer_Product);
MK_FUN_2(SYMMETRIC_TRANSPOSE_TIMES,Symmetric_Transpose_Times);
MK_FUN_1(TRANSPOSE_TIMES_SELF,Transpose_Times_Self);
MK_FUN_2(SYMMETRIC_TENSOR_PRODUCT_12,Symmetric_Tensor_Product_12);
MK_FUN_3(SYMMETRIC_DOUBLE_CONTRACT_12,Symmetric_Double_Contract_12);
MK_FUN_3(DOUBLE_CONTRACT_12,Double_Contract_12);

typedef MAT_MAP_1<NEG_B<ARG<0> >,VEC_NEG_B> MAT_NEG_B;
typedef MAT_MAP_2<ADD_BB<ARG<0>,ARG<1> >,VEC_ADD_BB> MAT_ADD_BB;
typedef MAT_MAP_2<SUB_BB<ARG<0>,ARG<1> >,VEC_SUB_BB> MAT_SUB_BB;
typedef MAT_MAP_1<MUL_BS<ARG<0>,ARG<1> >,VEC_MUL_BS> MAT_MUL_BS;
typedef MAT_MAP_1<DIV_BS<ARG<0>,ARG<1> >,VEC_DIV_BS> MAT_DIV_BS;
typedef OUTER_MAP<OUTER_PRODUCT_EQ<ARG<0> >,VEC_OUTER_PRODUCT_REV> MAT_OUTER_PRODUCT;
typedef SYM_OUTER_MAP<SYMMETRIC_OUTER_PRODUCT<ARG<0>,ARG<1> >,VEC_MAP_2<ADD_BB<OUTER_PRODUCT<ARG<3>,ARG<0> >,OUTER_PRODUCT<ARG<2>,ARG<1> > > > > MAT_SYM_OUTER_PRODUCT;
typedef SYM_OUTER_MAP<SYMMETRIC_TRANSPOSE_TIMES<ARG<0>,ARG<1> >,VEC_MAP_2<ADD_BB<TRANS_MUL<ARG<3>,ARG<0> >,TRANS_MUL<ARG<2>,ARG<1> > > > > MAT_SYM_TRANSPOSE_TIMES;
typedef OUTER_MAP<TRANSPOSE_TIMES_SELF<ARG<0> >,VEC_MAP_1<TRANS_MUL<ARG<1>,ARG<0> > > > MAT_TRANSPOSE_TIMES_SELF;
typedef MAT_MAP_1<TENSOR_PRODUCT_0<ARG<0>,ARG<1> >,TENSOR_PRODUCT_0_VM_V> MAT_TENSOR_PRODUCT_0;
typedef SYM_OUTER_MAP<SYMMETRIC_TENSOR_PRODUCT_12<ARG<0>,ARG<1> >,VEC_MAP_2<ADD_BB<TENSOR_PRODUCT_1<ARG<0>,ARG<3> >,TENSOR_PRODUCT_2<ARG<2>,ARG<1> > > > > MAT_SYM_TENSOR_PRODUCT_12;
typedef MAT_MAP_1<CONTRACT_00<ARG<0>,ARG<1> >,CONTRACT_00_VT_M> MAT_CONTRACT_00;
typedef MAT_MAP_1<CONTRACT_01<ARG<0>,ARG<1> >,CONTRACT_01_VT_M> MAT_CONTRACT_01;
typedef MAT_MAP_1<CONTRACT_0<ARG<0>,ARG<1> >,CONTRACT_0_VT_V> MAT_CONTRACT_0;
typedef SYM_OUTER_MAP<SYMMETRIC_DOUBLE_CONTRACT_12<ARG<2>,ARG<0>,ARG<1> >,VEC_MAP_2<SUB_BB<DOUBLE_CONTRACT_12<ARG<4>,ARG<2>,ARG<1> >,DOUBLE_CONTRACT_12<ARG<4>,ARG<3>,ARG<0> > > > > MAT_SYM_DOUBLE_CONTRACT_12;
typedef MAT_MAP_2<CHOOSE<ARG<0>,ARG<1> >,VEC_CHOOSE> MAT_CHOOSE;
typedef MAT_MAP_1<CHOOSE_ZERO<ARG<0> >,VEC_CHOOSE_ZERO> MAT_CHOOSE_ZERO;

}
}
#endif
