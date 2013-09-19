//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COL_TEN
//##################################################################### 
#ifndef __COL_TEN__
#define __COL_TEN__

#include <Tools/Auto_Diff/COL_MAT.h>
#include <Tools/Auto_Diff/PRIMITIVE_TENSORS.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV,class H> struct COL_TEN;

template<int nn,int m,int a,int b> struct CT_FLAGS
{
    enum WA{n=nn,mask=m,mask_a=a,mask_b=b};
    STATIC_ASSERT(n>=0 && mask<(1<<n) && mask_a<(1<<n) && mask_b<(1<<n));

    CT_FLAGS<n,m,a,b> operator+() const {return CT_FLAGS<n,m,a,b>();}
    CT_FLAGS<n,m,a,b> operator-() const {return CT_FLAGS<n,m,a,b>();}

    template<int m2,int a2,int b2>
    CT_FLAGS<n,(~m&~a&~b&m2)|(~m2&~a2&~b2&m)|(m&m2&~(a^a2)&~(b^b2)),(~m&~a&~b&a2)|(~m2&~a2&~b2&a)|(a&a2&~(m^m2)&~(b^b2)),b|b2|(m&m2&~b&~b2&(a^a2))>
    operator+(CT_FLAGS<n,m2,a2,b2> x) const
    {return CT_FLAGS<n,m&m2&~(a^a2)&~(b^b2),a&a2&~(m^m2)&~(b^b2),b|b2|(((a^a2)|(m^m2))&(m|a)&(m2|a2))>();}

    template<int m2,int a2,int b2>
    CT_FLAGS<n,(~m&~a&~b&m2)|(~m2&~a2&~b2&m)|(m&m2&~(a^a2)&~(b^b2)),(~m&~a&~b&a2)|(~m2&~a2&~b2&a)|(a&a2&~(m^m2)&~(b^b2)),b|b2|(m&m2&~b&~b2&(a^a2))>
    operator-(CT_FLAGS<n,m2,a2,b2> x) const
    {return CT_FLAGS<n,m&m2&~(a^a2)&~(b^b2),a&a2&~(m^m2)&~(b^b2),b|b2|(((a^a2)|(m^m2))&(m|a)&(m2|a2))>();}

    CT_FLAGS<n,m,a,b> operator*(double s) const {return CT_FLAGS<n,m,a,b>();}
    CT_FLAGS<n,m,a,b> operator/(double s) const {return CT_FLAGS<n,m,a,b>();}
    CT_FLAGS<n-1,(m>>1),(a>>1),(b>>1)> static Base() {return CT_FLAGS<n-1,(m>>1),(a>>1),(b>>1)>();}
    T_FLAGS<m&1,a&1,b&1> static First();
    template<int m2,int a2,int b2> CT_FLAGS<n+1,(m<<1)|m2,(a<<1)|a2,(b<<1)|b2> static Append(T_FLAGS<m2,a2,b2>);
};
typedef CT_FLAGS<0,0,0,0> CT_EMPTY;

#define MK_CT(TV,FLAGS) COL_TEN<TV,decltype(FLAGS)>

template<class TV>
struct COL_TEN<TV,CT_EMPTY>
{
    typedef CT_EMPTY FLAGS;
    typedef typename TV::SCALAR T;

    COL_TEN(){}

    void Neg(const COL_TEN& u) {}

    void Scale(const COL_TEN& u,T a) {}

    void Add(const COL_TEN& u,const COL_TEN& v) {}

    void Sub(const COL_TEN& u,const COL_TEN& v) {}

    template<class TV2>
    void Tensor_Product_0(const COL_MAT<TV,CM_EMPTY>& m,const TV2& v) {}

    template<class T_MAT>
    void Tensor_Product_1(const T_MAT& m,const COL_VEC<TV,CV_EMPTY>& v) {}

    template<class TV2>
    void Transpose_Tensor_Product_2(const COL_MAT<TV,CM_EMPTY>& m,const TV2& v) {}

    void Fill_From(const COL_TEN& u) {}
};

// Middle index is the long index
template<class TV,class R>
struct COL_TEN:public MK_CT(TV,R::Base())
{
    typedef R FLAGS;
    typedef typename TV::SCALAR T;

    typedef MK_TEN(TV,R::First()) T_DATA;

    T_DATA x;
    typedef MK_CT(TV,R::Base()) BASE;
    const BASE& Base() const {return static_cast<const BASE&>(*this);}
    BASE& Base() {return static_cast<BASE&>(*this);}
    typedef COL_TEN<TV,CT_FLAGS<R::n,0,0,0> > T_EMPTY;

    COL_TEN(){}

    template<class R2>
    void Neg(const COL_TEN<TV,R2>& u)
    {x=-u.x;Base().Neg(u.Base());}

    template<class R2>
    void Scale(const COL_TEN<TV,R2>& u,T a)
    {x=u.x*a;Base().Scale(u.Base(),a);}

    template<class R2,class R3>
    void Add(const COL_TEN<TV,R2>& u,const COL_TEN<TV,R3>& v)
    {x=u.x+v.x;Base().Add(u.Base(),v.Base());}

    template<class R2,class R3>
    void Sub(const COL_TEN<TV,R2>& u,const COL_TEN<TV,R3>& v)
    {x=u.x-v.x;Base().Sub(u.Base(),v.Base());}

    template<class TV2,class H>
    void Tensor_Product_0(const COL_MAT<TV,H>& m,const TV2& v)
    {x=::PhysBAM::HETERO_DIFF::Tensor_Product_0(m.x,v);Base().Tensor_Product_0(m.Base(),v);}

    template<class G,class T_MAT>
    void Tensor_Product_1(const T_MAT& m,const COL_VEC<TV,G>& v)
    {x=::PhysBAM::HETERO_DIFF::Tensor_Product_1(m,v.x);Base().Tensor_Product_1(m,v.Base());}

    template<class TV2,class H>
    void Transpose_Tensor_Product_2(const COL_MAT<TV,H>& m,const TV2& v)
    {x=::PhysBAM::HETERO_DIFF::Tensor_Product_2(m.x.Transposed(),v);Base().Transpose_Tensor_Product_2(m.Base(),v);}

    template<class R2>
    void Fill_From(const COL_TEN<TV,R2>& u)
    {Fill_From_Helper(x,u.x);Base().Fill_From(u.Base());}
};

template<class TV,class TV2> void
Contract_0(COL_MAT<TV,CM_EMPTY>& cm,const COL_TEN<TV,CT_EMPTY>& ct,const TV2& v) {}

template<class TV,class H,class R,class TV2> void
Contract_0(COL_MAT<TV,H>& cm,const COL_TEN<TV,R>& ct,const TV2& v)
{
    cm.x=Contract_0(ct.x,v);
    Contract_0(cm.Base(),ct.Base(),v);
}

template<class R> struct C0_CT_V;

template<> struct C0_CT_V<CT_EMPTY> {typedef CM_EMPTY TYPE;};

template<class R>
struct C0_CT_V
{
    typedef typename C0_CT_V<decltype(R::Base())>::TYPE B;
    typedef decltype(C0(R::First(),V_FLAGS<1>())) F;
    typedef decltype(B().Append(F())) TYPE;
};

template<int n,int m,int a,int b> CM_FLAGS<n,0,0> C0(CT_FLAGS<n,m,a,b>,V_FLAGS<0>);
template<int n,int m,int a,int b> typename C0_CT_V<CT_FLAGS<n,m,a,b> >::TYPE C0(CT_FLAGS<n,m,a,b>,V_FLAGS<1>);

template<class TV,class T_MAT> void
Contract_0(COL_TEN<TV,CT_EMPTY>& cm,const COL_TEN<TV,CT_EMPTY>& ct,const T_MAT& v) {}

template<class TV,class R2,class R,class T_MAT> void
Contract_0(COL_TEN<TV,R2>& ct2,const COL_TEN<TV,R>& ct,const T_MAT& m)
{
    ct2.x=Contract_0(ct.x,m);
    Contract_0(ct2.Base(),ct.Base(),m);
}

template<class R,class H> struct C0_CT_M;

template<class H> struct C0_CT_M<CT_EMPTY,H> {typedef CT_EMPTY TYPE;};

template<class R,class H>
struct C0_CT_M
{
    typedef typename C0_CT_M<decltype(R::Base()),H>::TYPE B;
    typedef decltype(C0(R::First(),H())) F;
    typedef decltype(B().Append(F())) TYPE;
};

template<int n,int m,int a,int b,int m2,int s2> typename C0_CT_M<CT_FLAGS<n,m,a,b>,M_FLAGS<m2,s2> >::TYPE C0(CT_FLAGS<n,m,a,b>,M_FLAGS<m2,s2>);

template<class TV,class T_TENS>
void Contract_1_Transpose(COL_TEN<TV,CT_EMPTY>& ct,const T_TENS& t,const COL_MAT<TV,CM_EMPTY>& cm) {}

template<class TV,class R,class H,class T_TENS>
void Contract_1_Transpose(COL_TEN<TV,R>& ct,const T_TENS& t,const COL_MAT<TV,H>& cm)
{
    ct.x=Contract_1(t,cm.x.Transposed());
    Contract_1_Transpose(ct.Base(),t,cm.Base());
}

template<class R,class H> struct C1t_T_CM;

template<class R> struct C1t_T_CM<R,CM_EMPTY> {typedef CT_EMPTY TYPE;};

template<class R,class H>
struct C1t_T_CM
{
    typedef typename C1t_T_CM<R,decltype(H::Base())>::TYPE B;
    typedef decltype(C1(R(),H::First())) F;
    typedef decltype(B().Append(F())) TYPE;
};

template<class R,class H> typename C1t_T_CM<R,H>::TYPE C1t(R,H);

}}
#endif
