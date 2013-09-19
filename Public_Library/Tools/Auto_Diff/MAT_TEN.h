//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MAT_TEN
//##################################################################### 
#ifndef __MAT_TEN__
#define __MAT_TEN__

#include <Tools/Auto_Diff/COL_TEN.h>
#include <Tools/Auto_Diff/MAT_MAT.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV,class H> struct MAT_TEN;

template<int nn,int m,int a,int b> struct MT_FLAGS
{
    enum WA{n=nn,mask=m,mask_a=a,mask_b=b};
    STATIC_ASSERT(mask<(1<<(n*(n+1)/2)) && mask_a<(1<<(n*(n+1)/2)) && mask_b<(1<<(n*(n+1)/2)));

    MT_FLAGS<n,m,a,b> operator+() const;
    MT_FLAGS<n,m,a,b> operator-() const;

    template<int m2,int a2,int b2>
    MT_FLAGS<n,(~m&~a&~b&m2)|(~m2&~a2&~b2&m)|(m&m2&~(a^a2)&~(b^b2)),(~m&~a&~b&a2)|(~m2&~a2&~b2&a)|(a&a2&~(m^m2)&~(b^b2)),b|b2|(m&m2&~b&~b2&(a^a2))>
    operator+(MT_FLAGS<n,m2,a2,b2> x) const;

    template<int m2,int a2,int b2>
    MT_FLAGS<n,(~m&~a&~b&m2)|(~m2&~a2&~b2&m)|(m&m2&~(a^a2)&~(b^b2)),(~m&~a&~b&a2)|(~m2&~a2&~b2&a)|(a&a2&~(m^m2)&~(b^b2)),b|b2|(m&m2&~b&~b2&(a^a2))>
    operator-(MT_FLAGS<n,m2,a2,b2> x) const;

    MT_FLAGS<n,m,a,b> operator*(double s) const;
    MT_FLAGS<n,m,a,b> operator/(double s) const;
    MT_FLAGS<n-1,(mask>>n),(mask_a>>n),(mask_b>>n)> static Base();
    CT_FLAGS<n-1,((mask&((1<<n)-1))>>1),((mask_a&((1<<n)-1))>>1),((mask_b&((1<<n)-1))>>1)> static Col();
    T_FLAGS<mask&1,mask_a&1,mask_b&1> static First();

    template<int m2,int a2,int b2,int m3,int a3,int b3>
    MT_FLAGS<n+1,m2|(m3<<1)|(m<<(n+1)),a2|(a3<<1)|(a<<(n+1)),b2|(b3<<1)|(b<<(n+1))>
    Append(T_FLAGS<m2,a2,b2>,CT_FLAGS<n,m3,a3,b3>);
};

typedef MT_FLAGS<0,0,0,0> MT_EMPTY;

template<int n,int m2,int m,int s> MT_FLAGS<n,s*m2,0,(m&~s)*m2> Tp0(V_FLAGS<m2>,MM_FLAGS<n,m,s>);

template<int n,int m,int m2,int s2>
CT_FLAGS<n,m&-s2,0,m&(-m2|-s2)> Tp1(CV_FLAGS<n,m>,M_FLAGS<m2,s2>);

template<int m,int n,int m2,int s2>
CT_FLAGS<n,-m&s2,-m&s2,-m&m2&~s2> Tp2(V_FLAGS<m>,CM_FLAGS<n,m2,s2>);

template<class G,class H> struct H_SYM_TP_12;
template<> struct H_SYM_TP_12<CV_EMPTY,CM_EMPTY> {typedef MT_EMPTY TYPE;};
template<class G,class H> struct H_SYM_TP_12
{
    typedef decltype(Tp1(G::Base(),H::First())+Tp2(G::First(),H::Base())) COL;
    typedef typename H_SYM_TP_12<decltype(G::Base()),decltype(H::Base())>::TYPE B;
    typedef decltype(Stp12(G::First(),H::First())) M;
    typedef decltype(B().Append(M(),COL())) TYPE;
};
template<int n,int m,int m2,int s2> typename H_SYM_TP_12<CV_FLAGS<n,m>,CM_FLAGS<n,m2,s2> >::TYPE
Stp12(CV_FLAGS<n,m>,CM_FLAGS<n,m2,s2>);

template<int n,int m,int a,int b> MM_FLAGS<n,0,0> C0(MT_FLAGS<n,m,a,b>,V_FLAGS<0>);
template<int n,int m,int a,int b> MM_FLAGS<n,m|a|b,m&~a&~b> C0(MT_FLAGS<n,m,a,b>,V_FLAGS<1>);

#define MK_MT(TV,FLAGS) MAT_TEN<TV,decltype(FLAGS)>

template<class TV>
struct MAT_TEN<TV,MT_EMPTY>
{
    typedef MT_EMPTY FLAGS;
    typedef typename TV::SCALAR T;

    MATRIX<T,TV::m> operator()(int i,int j) const {return MATRIX<T,TV::m>();}

    SYMMETRIC_MATRIX<T,TV::m> operator()(int i) const {return SYMMETRIC_MATRIX<T,TV::m>();}

    void Neg(const MAT_TEN& u) {}

    void Scale(const MAT_TEN& u,T a) {}

    void Add(const MAT_TEN& u,const MAT_TEN& v) {}

    void Sub(const MAT_TEN& u,const MAT_TEN& v) {}

    void Tensor_Product_0(const MAT_MAT<TV,MM_EMPTY>& m,const TV& v) {}

    void Fill_From(const MAT_TEN& u) {}
};

template<class TV,class R>
struct MAT_TEN:public MK_MT(TV,R::Base())
{
    typedef R FLAGS;
    typedef typename TV::SCALAR T;

    typedef MK_CT(TV,R::Col()) T_COL;
    typedef MK_SYM_TEN(TV,R::First()) T_DATA;
    T_DATA x;
    T_COL c;

    typedef MK_MT(TV,R::Base()) BASE;
    const BASE& Base() const {return static_cast<const BASE&>(*this);}
    BASE& Base() {return static_cast<BASE&>(*this);}

    template<class R2>
    void Neg(const MAT_TEN<TV,R2>& u)
    {x=-u.x;c.Neg(u.c);Base().Neg(u.Base());}

    template<class R2>
    void Scale(const MAT_TEN<TV,R2>& u,T a)
    {x=u.x*a;c.Scale(u.c,a);Base().Scale(u.Base(),a);}

    template<class R2,class R3>
    void Add(const MAT_TEN<TV,R2>& u,const MAT_TEN<TV,R3>& v)
    {x=u.x+v.x;c.Add(u.c,v.c);Base().Add(u.Base(),v.Base());}

    template<class R2,class R3>
    void Sub(const MAT_TEN<TV,R2>& u,const MAT_TEN<TV,R3>& v)
    {x=u.x-v.x;c.Sub(u.c,v.c);Base().Sub(u.Base(),v.Base());}

    template<class H>
    void Tensor_Product_0(const MAT_MAT<TV,H>& m,const TV& v)
    {x=::PhysBAM::HETERO_DIFF::Tensor_Product_0(m.x,v);c.Tensor_Product_0(m.c,v);Base().Tensor_Product_0(m.Base(),v);}

    template<class R2>
    void Fill_From(const MAT_TEN<TV,R2>& u)
    {Fill_From_Helper(x,u.x);c.Fill_From(u.c);Base().Fill_From(u.Base());}
};

template<class TV> void
Contract_0(MAT_MAT<TV,MM_EMPTY>& mm,const MAT_TEN<TV,MT_EMPTY>& mt,const TV& v) {}

template<class TV,class H,class R> void
Contract_0(MAT_MAT<TV,H>& mm,const MAT_TEN<TV,R>& mt,const TV& v)
{
    mm.x=Contract_0(mt.x,v);
    Contract_0(mm.c,mt.c,v);
    Contract_0(mm.Base(),mt.Base(),v);
}

template<class TV> void
Symmetric_Double_Contract_12_With_Transposes(MAT_TEN<TV,MT_EMPTY>& mt,const PERM_TENSOR<TV>&t,const COL_MAT<TV,CM_EMPTY>& cm1,const COL_MAT<TV,CM_EMPTY>& cm2) {}

template<class TV,class H1,class H2,class R> void
Symmetric_Double_Contract_12_With_Transposes(MAT_TEN<TV,R>& mt,const PERM_TENSOR<TV>&t,const COL_MAT<TV,H1>& cm1,const COL_MAT<TV,H2>& cm2)
{
    auto t2_cm1=Contract_2(t,cm1.x.Transposed());
    auto t2_cm2=Contract_2(t,cm2.x.Transposed());
    mt.x=Symmetric_Double_Contract_12(t,cm1.x.Transposed(),cm2.x.Transposed());

    typedef typename GET_FLAGS<decltype(t2_cm1)>::TYPE CM1;
    typedef typename GET_FLAGS<decltype(t2_cm2)>::TYPE CM2;

    MK_CT(TV,C1t(CM2(),H1::Base())) c0;
    MK_CT(TV,C1t(CM1(),H2::Base())) c1;
    Contract_1_Transpose(c0,t2_cm2,cm1.Base());
    Contract_1_Transpose(c1,t2_cm1,cm2.Base());
    mt.c.Sub(c0,c1);
    Symmetric_Double_Contract_12_With_Transposes(mt.Base(),t,cm1.Base(),cm2.Base());
}

template<class TV> void
Contract_0(MAT_TEN<TV,MT_EMPTY>& mt2,const MAT_TEN<TV,MT_EMPTY>& mt,const MATRIX<typename TV::SCALAR,TV::m>& m) {}

template<class TV,class R,class R2> void
Contract_0(MAT_TEN<TV,R2>& mt2,const MAT_TEN<TV,R>& mt,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    mt2.x=Contract_0(mt.x,m);
    Contract_0(mt2.c,mt.c,m);
    Contract_0(mt2.Base(),mt.Base(),m);
}

template<class R> struct C0_MT_M;

template<> struct C0_MT_M<MT_EMPTY> {typedef MT_EMPTY TYPE;};

template<class R>
struct C0_MT_M
{
    typedef typename C0_MT_M<decltype(R::Base())>::TYPE B;
    typedef decltype(C0(R::First(),M_FLAGS<1,0>())) F;
    typedef decltype(C0(R::Col(),M_FLAGS<1,0>())) C;
    typedef decltype(B().Append(F(),C())) TYPE;
};

template<int n,int m,int a,int b> typename C0_MT_M<MT_FLAGS<n,m,a,b> >::TYPE C0(MT_FLAGS<n,m,a,b>,M_FLAGS<1,0>);


}
}

#endif
