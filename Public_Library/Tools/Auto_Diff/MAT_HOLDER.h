//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MAT_HOLDER
//##################################################################### 
#ifndef __MAT_HOLDER__
#define __MAT_HOLDER__

#include <Tools/Auto_Diff/VEC_HOLDER.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

struct MAT_END {};

template<class OBJ_IN,class COL_IN,class BASE_IN>
struct MAT_HOLDER
{
    typedef typename OBJ_IN::SCALAR T;
    typedef OBJ_IN OBJ;
    typedef COL_IN COL;
    typedef BASE_IN BASE;
    OBJ x;
    COL y;
    BASE z;
};

inline void Fill_From(MAT_END& o,const MAT_END& m) {}
template<class OBJ,class COL,class BASE,class OBJ2,class COL2,class BASE2>
void Fill_From(MAT_HOLDER<OBJ2,COL2,BASE2>& o,const MAT_HOLDER<OBJ,COL,BASE>& m)
{
    Fill_From(o.x,m.x);
    Fill_From(o.y,m.y);
    Fill_From(o.z,m.z);
}

template<class OUT,class OBJ,class COL,class BASE>
void Get(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& m,int i,int j)
{
    if(i<j){Get_Helper(o,m,j,i);o=o.Transposed();}
    else Get_Helper(o,m,i,j);
}
template<class OUT> void Get_Helper(OUT& o,const MAT_END& v,int i,int j) {PHYSBAM_FATAL_ERROR();}
template<class OUT,class OBJ,class COL,class BASE>
void Get_Helper(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& m,int i,int j)
{
    if(j>0) return Get_Helper(o,m.z,i-1,j-1);
    if(i>0) return Get(o,m.y,i-1);
    Fill_From(o,m.x);
}

template<class OUT> void Get_Diag(OUT& o,const MAT_END& v,int i) {PHYSBAM_FATAL_ERROR();}
template<class OUT,class OBJ,class COL,class BASE>
void Get_Diag(OUT& o,const MAT_HOLDER<OBJ,COL,BASE>& m,int i)
{
    if(i>0) return Get_Diag(o,m.z,i-1);
    Fill_From(o,m.x);
}

template<class TV,int n> struct EMPTY_MAT
{typedef MAT_HOLDER<ZERO_MAT<TV>,typename EMPTY_VEC_MAT<TV,n-1>::TYPE,typename EMPTY_MAT<TV,n-1>::TYPE> TYPE;};
template<class TV> struct EMPTY_MAT<TV,0> {typedef MAT_END TYPE;};

template<class TV,int n> struct EMPTY_MAT_TEN
{typedef MAT_HOLDER<ZERO_TENSOR<TV>,typename EMPTY_VEC_TEN<TV,n-1>::TYPE,typename EMPTY_MAT_TEN<TV,n-1>::TYPE> TYPE;};
template<class TV> struct EMPTY_MAT_TEN<TV,0> {typedef MAT_END TYPE;};

template<class OP,class VEC_OP,class MAT> struct TYPE_MAT_MAP_1;

template<class OP,class VEC_OP>
struct TYPE_MAT_MAP_1<OP,VEC_OP,MAT_END>
{template<class ...Args> static MAT_END Type(Args&&...);};

template<class OP,class VEC_OP,class OBJ,class COL,class BASE>
struct TYPE_MAT_MAP_1<OP,VEC_OP,MAT_HOLDER<OBJ,COL,BASE> >
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(VEC_OP::Type(COL(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(TYPE_MAT_MAP_1<OP,VEC_OP,BASE>::Type(typename REMOVE_REFERENCE<Args>::TYPE()...))> Type(Args&&...);};

template<class OP,class VEC_OP> struct MAT_MAP_1
{
    template<class OBJ,class COL,class BASE,class ...Args>
    static decltype(TYPE_MAT_MAP_1<OP,VEC_OP,MAT_HOLDER<OBJ,COL,BASE> >::Type(typename REMOVE_REFERENCE<Args>::TYPE()...))
    Type(const MAT_HOLDER<OBJ,COL,BASE>&,Args&&...);
    template<class ...Args> static MAT_END Type(MAT_END,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const MAT_END& in,Args&&... args) const {}
    template<class OBJ,class COL,class BASE,class ...Args> void operator()(decltype(Type(MAT_HOLDER<OBJ,COL,BASE>(),typename REMOVE_REFERENCE<Args>::TYPE()...))& out,const MAT_HOLDER<OBJ,COL,BASE>& in,Args&&... args) const
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
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ0(),OBJ1(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(VEC_OP::Type(COL0(),COL1(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(TYPE_MAT_MAP_2<OP,VEC_OP,BASE0,BASE1>::Type(typename REMOVE_REFERENCE<Args>::TYPE()...))>Type(Args&&...);};

template<class OP,class VEC_OP> struct MAT_MAP_2
{
    template<class OBJ0,class COL0,class BASE0,class OBJ1,class COL1,class BASE1,class ...Args>
    static decltype(TYPE_MAT_MAP_2<OP,VEC_OP,MAT_HOLDER<OBJ0,COL0,BASE0>,MAT_HOLDER<OBJ1,COL1,BASE1> >::Type(typename REMOVE_REFERENCE<Args>::TYPE()...))
    Type(const MAT_HOLDER<OBJ0,COL0,BASE0>&,const MAT_HOLDER<OBJ1,COL1,BASE1>&,Args&&...);
    template<class ...Args> static MAT_END Type(MAT_END,MAT_END,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const MAT_END& in0,const MAT_END& in1,Args&&... args) const {}
    template<class OBJ0,class COL0,class BASE0,class OBJ1,class COL1,class BASE1,class ...Args> void
    operator()(decltype(Type(MAT_HOLDER<OBJ0,COL0,BASE0>(),MAT_HOLDER<OBJ1,COL1,BASE1>(),typename REMOVE_REFERENCE<Args>::TYPE()...))& out,const MAT_HOLDER<OBJ0,COL0,BASE0>& in0,const MAT_HOLDER<OBJ1,COL1,BASE1>& in1,Args&&... args) const
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
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ0(),OBJ1(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(VEC_OP::Type(BASE0(),BASE1(),OBJ0(),OBJ1(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(TYPE_SYM_OUTER_MAP<OP,VEC_OP,BASE0,BASE1>::Type(typename REMOVE_REFERENCE<Args>::TYPE()...))> Type(Args&&...);};

template<class OP,class VEC_OP> struct SYM_OUTER_MAP
{
    template<class ...Args> static MAT_END Type(VEC_END,VEC_END,Args&&... args);
    template<class OBJ0,class BASE0,class OBJ1,class BASE1,class ...Args> static
    decltype(TYPE_SYM_OUTER_MAP<OP,VEC_OP,VEC_HOLDER<OBJ0,BASE0>,VEC_HOLDER<OBJ1,BASE1> >::Type(typename REMOVE_REFERENCE<Args>::TYPE()...)) Type(const VEC_HOLDER<OBJ0,BASE0>&,const VEC_HOLDER<OBJ1,BASE1>&,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const VEC_END& in0,const VEC_END& in1,Args&&... args) const {}
    template<class OBJ0,class BASE0,class OBJ1,class BASE1,class ...Args> void
    operator()(decltype(Type(VEC_HOLDER<OBJ0,BASE0>(),VEC_HOLDER<OBJ1,BASE1>(),typename REMOVE_REFERENCE<Args>::TYPE()...))& out,const VEC_HOLDER<OBJ0,BASE0>& in0,const VEC_HOLDER<OBJ1,BASE1>& in1,Args&&... args) const
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
{template<class ...Args> static MAT_HOLDER<decltype(OP()(OBJ(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(VEC_OP::Type(BASE(),OBJ(),typename REMOVE_REFERENCE<Args>::TYPE()...)),decltype(TYPE_OUTER_MAP<OP,VEC_OP,BASE>::Type(typename REMOVE_REFERENCE<Args>::TYPE()...))> Type(Args&&...);};

template<class OP,class VEC_OP> struct OUTER_MAP
{
    template<class ...Args> static MAT_END Type(VEC_END,VEC_END,Args&&... args);
    template<class OBJ,class BASE,class ...Args> static
    decltype(TYPE_OUTER_MAP<OP,VEC_OP,VEC_HOLDER<OBJ,BASE> >::Type(typename REMOVE_REFERENCE<Args>::TYPE()...)) Type(const VEC_HOLDER<OBJ,BASE>&,Args&&...);

    template<class ...Args> void operator()(MAT_END& out,const VEC_END& in,Args&&... args) const {}
    template<class OBJ,class BASE,class ...Args> void
    operator()(decltype(Type(VEC_HOLDER<OBJ,BASE>(),typename REMOVE_REFERENCE<Args>::TYPE()...))& out,const VEC_HOLDER<OBJ,BASE>& in,Args&&... args) const
    {
        out.x=OP()(in.x,args...);
        VEC_OP()(out.y,in.z,in.x,args...);
        (*this)(out.z,in.z,args...);
    };
};

MK_FUN_1(OUTER_PRODUCT_HELPER_EQ,Outer_Product_Helper);
MK_FUN_2(SYMMETRIC_OUTER_PRODUCT_HELPER,Symmetric_Outer_Product_Helper);
MK_FUN_2(SYMMETRIC_TRANSPOSE_TIMES,Symmetric_Transpose_Times);
MK_FUN_1(TRANSPOSE_TIMES_SELF,Transpose_Times_Self);
MK_FUN_2(SYMMETRIC_TENSOR_PRODUCT_12,Symmetric_Tensor_Product_12);
MK_FUN_3(SYMMETRIC_DOUBLE_CONTRACT_12,Symmetric_Double_Contract_12);

typedef MAT_MAP_1<NEG<ARG<0> >,VEC_NEG> MAT_NEG;
typedef MAT_MAP_2<ADD<ARG<0>,ARG<1> >,VEC_ADD> MAT_ADD;
typedef MAT_MAP_2<SUB<ARG<0>,ARG<1> >,VEC_SUB> MAT_SUB;
typedef MAT_MAP_1<MUL<ARG<0>,ARG<1> >,VEC_SCALE> MAT_SCALE;
typedef MAT_MAP_1<DIV<ARG<0>,ARG<1> >,VEC_SCALE_DIV> MAT_SCALE_DIV;
typedef OUTER_MAP<OUTER_PRODUCT_HELPER_EQ<ARG<0> >,VEC_OUTER_PRODUCT> MAT_OUTER_PRODUCT;
typedef SYM_OUTER_MAP<SYMMETRIC_OUTER_PRODUCT_HELPER<ARG<0>,ARG<1> >,VEC_MAP_2<ADD<OUTER_PRODUCT_HELPER<ARG<0>,ARG<3> >,OUTER_PRODUCT_HELPER<ARG<1>,ARG<2> > > > > MAT_SYM_OUTER_PRODUCT;
typedef SYM_OUTER_MAP<SYMMETRIC_TRANSPOSE_TIMES<ARG<0>,ARG<1> >,VEC_MAP_2<ADD<MUL<TRANSPOSE<ARG<0> >,ARG<3> >,MUL<TRANSPOSE<ARG<1> >,ARG<2> > > > > MAT_SYM_TRANSPOSE_TIMES;
typedef OUTER_MAP<TRANSPOSE_TIMES_SELF<ARG<0> >,VEC_MAP_1<MUL<TRANSPOSE<ARG<0> >,ARG<1> > > > MAT_TRANSPOSE_TIMES_SELF;
typedef MAT_MAP_1<TENSOR_PRODUCT_0<ARG<0>,ARG<1> >,TENSOR_PRODUCT_0_VM_V> MAT_TENSOR_PRODUCT_0;
typedef SYM_OUTER_MAP<SYMMETRIC_TENSOR_PRODUCT_12<ARG<0>,ARG<1> >,VEC_MAP_2<ADD<TENSOR_PRODUCT_2<ARG<0>,ARG<3> >,TENSOR_PRODUCT_1<ARG<2>,ARG<1> > > > > MAT_SYM_TENSOR_PRODUCT_12;
typedef MAT_MAP_1<CONTRACT_0<ARG<0>,ARG<1> >,CONTRACT_0_VT_X> MAT_CONTRACT_0;
typedef SYM_OUTER_MAP<SYMMETRIC_DOUBLE_CONTRACT_12<ARG<2>,ARG<0>,ARG<1> >,VEC_MAP_2<SUB<CONTRACT_1<CONTRACT_2<ARG<4>,ARG<3> >,ARG<0> >,CONTRACT_1<CONTRACT_2<ARG<4>,ARG<2> >,ARG<1> > > > > MAT_SYM_DOUBLE_CONTRACT_12;
typedef MAT_MAP_2<CHOOSE<ARG<0>,ARG<1> >,VEC_CHOOSE> MAT_CHOOSE;

}
}
#endif
