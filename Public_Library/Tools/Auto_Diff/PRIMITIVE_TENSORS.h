//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRIMITIVE_TENSORS
//##################################################################### 
#ifndef __PRIMITIVE_TENSORS__
#define __PRIMITIVE_TENSORS__

#include <Tools/Auto_Diff/PRIMITIVE_MATRICES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<int m,int a,int b> struct T_FLAGS;

// symmetric T_ijk=T_ikj
template<class TV>
struct SYMMETRIC_TENSOR
{
    typedef T_FLAGS<0,0,1> FLAGS;
    typedef typename TV::SCALAR T;
    VECTOR<SYMMETRIC_MATRIX<T,TV::m>,TV::m> x;

    SYMMETRIC_TENSOR operator-() const
    {SYMMETRIC_TENSOR t;for(int i=0;i<TV::m;i++) t.x(i)=-x(i);return t;}

    SYMMETRIC_TENSOR operator*(T a) const
    {SYMMETRIC_TENSOR t;for(int i=0;i<TV::m;i++) t.x(i)=x(i)*a;return t;}

    SYMMETRIC_TENSOR operator/(T a) const
    {SYMMETRIC_TENSOR t;for(int i=0;i<TV::m;i++) t.x(i)=x(i)/a;return t;}
};

// general T_ijk
template<class TV>
struct TENSOR
{
    typedef T_FLAGS<0,0,1> FLAGS;
    typedef typename TV::SCALAR T;
    VECTOR<MATRIX<T,TV::m>,TV::m> x;

    TENSOR(){}

    TENSOR(const SYMMETRIC_TENSOR<TV>& t)
    {for(int i=0;i<TV::m;i++) x(i)=t.x(i);}

    TENSOR operator-() const
    {TENSOR t;for(int i=0;i<TV::m;i++) t.x(i)=-x(i);return t;}

    TENSOR operator*(T a) const
    {TENSOR t;for(int i=0;i<TV::m;i++) t.x(i)=x(i)*a;return t;}

    TENSOR operator/(T a) const
    {TENSOR t;for(int i=0;i<TV::m;i++) t.x(i)=x(i)/a;return t;}

    SYMMETRIC_TENSOR<TV> Twice_Symmetric_Part_12() const
    {SYMMETRIC_TENSOR<TV> t;for(int i=0;i<TV::m;i++) t.x(i)=x(i).Twice_Symmetric_Part();return t;}
};

// T_ijk=0
template<class TV>
struct ZERO_TENSOR
{
    typedef T_FLAGS<0,0,0> FLAGS;
    typedef typename TV::SCALAR T;
    ZERO_TENSOR operator-() const
    {return *this;}

    ZERO_TENSOR operator*(T a) const
    {return *this;}

    ZERO_TENSOR operator/(T a) const
    {return *this;}
};

// T_ijk=v_i delta_jk
template<class TV>
struct VEC_ID_TENSOR_0
{
    typedef T_FLAGS<1,0,0> FLAGS;
    TV v;
    explicit VEC_ID_TENSOR_0(TV v=TV()): v(v) {}

    typedef typename TV::SCALAR T;
    VEC_ID_TENSOR_0 operator-() const
    {return VEC_ID_TENSOR_0(-v);}

    VEC_ID_TENSOR_0 operator*(T a) const
    {return VEC_ID_TENSOR_0(v*a);}

    VEC_ID_TENSOR_0 operator/(T a) const
    {return VEC_ID_TENSOR_0(v/a);}
};

// T_ijk=v_j*delta_ik
template<class TV>
struct VEC_ID_TENSOR_1
{
    typedef T_FLAGS<1,0,1> FLAGS;
    TV v;
    explicit VEC_ID_TENSOR_1(TV v=TV()): v(v) {}

    typedef typename TV::SCALAR T;
    VEC_ID_TENSOR_1 operator-() const
    {return VEC_ID_TENSOR_1(-v);}

    VEC_ID_TENSOR_1 operator*(T a) const
    {return VEC_ID_TENSOR_1(v*a);}

    VEC_ID_TENSOR_1 operator/(T a) const
    {return VEC_ID_TENSOR_1(v/a);}
};

// T_ijk=v_k*delta_ij
template<class TV>
struct VEC_ID_TENSOR_2
{
    typedef T_FLAGS<1,1,0> FLAGS;
    TV v;
    explicit VEC_ID_TENSOR_2(TV v=TV()): v(v) {}

    typedef typename TV::SCALAR T;
    VEC_ID_TENSOR_2 operator-() const
    {return VEC_ID_TENSOR_2(-v);}

    VEC_ID_TENSOR_2 operator*(T a) const
    {return VEC_ID_TENSOR_2(v*a);}

    VEC_ID_TENSOR_2 operator/(T a) const
    {return VEC_ID_TENSOR_2(v/a);}
};

// T_ijk=v_j*delta_ik+v_k*delta_ij
template<class TV>
struct VEC_ID_TENSOR_12
{
    typedef T_FLAGS<1,1,1> FLAGS;
    TV v;
    explicit VEC_ID_TENSOR_12(TV v=TV()): v(v) {}

    typedef typename TV::SCALAR T;
    VEC_ID_TENSOR_12 operator-() const
    {return VEC_ID_TENSOR_12(-v);}

    VEC_ID_TENSOR_12 operator*(T a) const
    {return VEC_ID_TENSOR_12(v*a);}

    VEC_ID_TENSOR_12 operator/(T a) const
    {return VEC_ID_TENSOR_12(v/a);}
};

// T_ijk=x*e_ijk
template<class TV>
struct PERM_TENSOR
{
    typedef T_FLAGS<0,1,1> FLAGS;
    STATIC_ASSERT(TV::m==3);

    typedef typename TV::SCALAR T;
    T x;

    explicit PERM_TENSOR(T x=T()): x(x) {}

    PERM_TENSOR operator-() const
    {return PERM_TENSOR(-x);}

    PERM_TENSOR operator*(T a) const
    {return PERM_TENSOR(x*a);}

    PERM_TENSOR operator/(T a) const
    {return PERM_TENSOR(x/a);}
};

template<int m,int a,int b> struct T_FLAGS
{
    enum WA{mask=m,mask_a=a,mask_b=b};

    T_FLAGS<m,a,b> operator+() const {return T_FLAGS<m,a,b>();}
    T_FLAGS<m,a,b> operator-() const {return T_FLAGS<m,a,b>();}

    template<int m2,int a2,int b2>
    T_FLAGS<m&m2&~(a^a2)&~(b^b2),a&a2&~(m^m2)&~(b^b2),b|b2|(((a^a2)|(m^m2))&(m|a)&(m2|a2))>
    operator+(T_FLAGS<m2,a2,b2> x) const
    {return T_FLAGS<m&m2&~(a^a2)&~(b^b2),a&a2&~(m^m2)&~(b^b2),b|b2|(((a^a2)|(m^m2))&(m|a)&(m2|a2))>();}

    template<int m2,int a2,int b2>
    T_FLAGS<m&m2&~(a^a2)&~(b^b2),a&a2&~(m^m2)&~(b^b2),b|b2|(((a^a2)|(m^m2))&(m|a)&(m2|a2))>
    operator-(T_FLAGS<m2,a2,b2> x) const
    {return T_FLAGS<m&m2&~(a^a2)&~(b^b2),a&a2&~(m^m2)&~(b^b2),b|b2|(((a^a2)|(m^m2))&(m|a)&(m2|a2))>();}

    T_FLAGS<m,a,b> operator*(double s) const {return T_FLAGS<m,a,b>();}
    T_FLAGS<m,a,b> operator/(double s) const {return T_FLAGS<m,a,b>();}
};

template<class TV,class R> struct CH_TEN;
template<class TV> struct CH_TEN<TV,T_FLAGS<0,0,0> > {typedef ZERO_TENSOR<TV> TYPE;};
template<class TV> struct CH_TEN<TV,T_FLAGS<0,0,1> > {typedef TENSOR<TV> TYPE;};
template<class TV> struct CH_TEN<TV,T_FLAGS<0,1,1> > {typedef PERM_TENSOR<TV> TYPE;};
template<class TV> struct CH_TEN<TV,T_FLAGS<1,0,0> > {typedef VEC_ID_TENSOR_0<TV> TYPE;};
template<class TV> struct CH_TEN<TV,T_FLAGS<1,0,1> > {typedef VEC_ID_TENSOR_1<TV> TYPE;};
template<class TV> struct CH_TEN<TV,T_FLAGS<1,1,0> > {typedef VEC_ID_TENSOR_2<TV> TYPE;};
template<class TV> struct CH_TEN<TV,T_FLAGS<1,1,1> > {typedef VEC_ID_TENSOR_12<TV> TYPE;};

template<class TV,class R> struct CH_SYM_TEN;
template<class TV> struct CH_SYM_TEN<TV,T_FLAGS<0,0,0> > {typedef ZERO_TENSOR<TV> TYPE;};
template<class TV> struct CH_SYM_TEN<TV,T_FLAGS<0,0,1> > {typedef SYMMETRIC_TENSOR<TV> TYPE;};
template<class TV> struct CH_SYM_TEN<TV,T_FLAGS<1,0,0> > {typedef VEC_ID_TENSOR_0<TV> TYPE;};
template<class TV> struct CH_SYM_TEN<TV,T_FLAGS<1,1,1> > {typedef VEC_ID_TENSOR_12<TV> TYPE;};

#define MK_TEN(TV,R) typename CH_TEN<TV,decltype(R)>::TYPE
#define MK_SYM_TEN(TV,R) typename CH_SYM_TEN<TV,decltype(R)>::TYPE

template<class TV,class TN>
ZERO_MAT<TV> Contract_0(const TN& t,ZERO_VEC<TV> z){return ZERO_MAT<TV>();}

template<class TV>
MATRIX<typename TV::SCALAR,TV::m> Contract_0(const TENSOR<TV>& t,const TV& v)
{
    MATRIX<typename TV::SCALAR,TV::m> m;
    for(int i=0;i<TV::m;i++) m+=t.x(i)*v(i);
    return m;
}

template<class TV>
SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Contract_0(const SYMMETRIC_TENSOR<TV>& t,const TV& v)
{
    SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> m;
    for(int i=0;i<TV::m;i++) m+=t.x(i)*v(i);
    return m;
}

template<class TV>
ZERO_MAT<TV> Contract_0(const ZERO_TENSOR<TV>& t,const TV& v) {return ZERO_MAT<TV>();}

template<class TV>
SCALE_MAT<TV> Contract_0(const VEC_ID_TENSOR_0<TV>& t,const TV& v) {return SCALE_MAT<TV>(t.v.Dot(v));}

template<class TV>
MATRIX<typename TV::SCALAR,TV::m> Contract_0(const VEC_ID_TENSOR_1<TV>& t,const TV& v) {return MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(t.v,v);}

template<class TV>
MATRIX<typename TV::SCALAR,TV::m> Contract_0(const VEC_ID_TENSOR_2<TV>& t,const TV& v) {return MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(v,t.v);}

template<class TV>
SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> Contract_0(const VEC_ID_TENSOR_12<TV>& t,const TV& v) {return SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(t.v,v);}

template<class TV>
MATRIX<typename TV::SCALAR,TV::m> Contract_0(const PERM_TENSOR<TV>& t,const TV& v) {return MATRIX<typename TV::SCALAR,TV::m>::Cross_Product_Matrix(-t.x*v);}

template<int m,int a,int b> M_FLAGS<0,0> C0(T_FLAGS<m,a,b>,V_FLAGS<0>);
M_FLAGS<1,0> C0(T_FLAGS<0,0,1>,V_FLAGS<1>);
M_FLAGS<0,0> C0(T_FLAGS<0,0,0>,V_FLAGS<1>);
M_FLAGS<1,1> C0(T_FLAGS<1,0,0>,V_FLAGS<1>);
M_FLAGS<1,0> C0(T_FLAGS<1,0,1>,V_FLAGS<1>);
M_FLAGS<1,0> C0(T_FLAGS<1,1,0>,V_FLAGS<1>);
M_FLAGS<1,0> C0(T_FLAGS<1,1,1>,V_FLAGS<1>);
M_FLAGS<1,0> C0(T_FLAGS<0,1,1>,V_FLAGS<1>);

template<class TV,class MAT>
ZERO_TENSOR<TV> Tensor_Product_0(const MAT& t,ZERO_VEC<TV> z){return ZERO_TENSOR<TV>();}

template<class TV>
TENSOR<TV> Tensor_Product_0(const MATRIX<typename TV::SCALAR,TV::m>& m,const TV& v)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=m*v(i);
    return t;
}

template<class TV>
SYMMETRIC_TENSOR<TV> Tensor_Product_0(const SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>& m,const TV& v)
{
    SYMMETRIC_TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=m*v(i);
    return t;
}

template<class TV>
ZERO_TENSOR<TV> Tensor_Product_0(const ZERO_MAT<TV>& m,const TV& v) {return ZERO_TENSOR<TV>();}

template<class TV>
VEC_ID_TENSOR_0<TV> Tensor_Product_0(ID_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_0<TV>(v);}

template<class TV>
VEC_ID_TENSOR_0<TV> Tensor_Product_0(SCALE_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_0<TV>(m.x*v);}


template<class TV,class MAT>
ZERO_TENSOR<TV> Tensor_Product_1(const MAT& t,ZERO_VEC<TV> z){return ZERO_TENSOR<TV>();}

template<class TV>
TENSOR<TV> Tensor_Product_1(const MATRIX<typename TV::SCALAR,TV::m>& m,const TV& v)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(v,m.Row(i));
    return t;
}

template<class TV>
ZERO_TENSOR<TV> Tensor_Product_1(const ZERO_MAT<TV>& m,const TV& v) {return ZERO_TENSOR<TV>();}

template<class TV>
VEC_ID_TENSOR_1<TV> Tensor_Product_1(ID_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_1<TV>(v);}

template<class TV>
VEC_ID_TENSOR_1<TV> Tensor_Product_1(SCALE_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_1<TV>(m.x*v);}


template<class TV,class MAT>
ZERO_TENSOR<TV> Tensor_Product_2(const MAT& t,ZERO_VEC<TV> z){return ZERO_TENSOR<TV>();}

template<class TV>
TENSOR<TV> Tensor_Product_2(const MATRIX<typename TV::SCALAR,TV::m>& m,const TV& v)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(m.Row(i),v);
    return t;
}

template<class TV>
ZERO_TENSOR<TV> Tensor_Product_2(const ZERO_MAT<TV>& m,const TV& v) {return ZERO_TENSOR<TV>();}

template<class TV>
VEC_ID_TENSOR_2<TV> Tensor_Product_2(ID_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_2<TV>(v);}

template<class TV>
VEC_ID_TENSOR_2<TV> Tensor_Product_2(SCALE_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_2<TV>(m.x*v);}


template<class TV,class MAT>
ZERO_TENSOR<TV> Symmetric_Tensor_Product_12(const MAT& t,ZERO_VEC<TV> z){return ZERO_TENSOR<TV>();}

template<class TV>
SYMMETRIC_TENSOR<TV> Symmetric_Tensor_Product_12(const MATRIX<typename TV::SCALAR,TV::m>& m,const TV& v)
{
    SYMMETRIC_TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(m.Row(i),v);
    return t;
}

template<class TV>
ZERO_TENSOR<TV> Symmetric_Tensor_Product_12(const ZERO_MAT<TV>& m,const TV& v) {return ZERO_TENSOR<TV>();}

template<class TV>
VEC_ID_TENSOR_12<TV> Symmetric_Tensor_Product_12(ID_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_12<TV>(v);}

template<class TV>
VEC_ID_TENSOR_12<TV> Symmetric_Tensor_Product_12(SCALE_MAT<TV> m,const TV& v) {return VEC_ID_TENSOR_12<TV>(m.x*v);}

template<int m2,int s2> T_FLAGS<0,0,0> Stp12(V_FLAGS<0>,M_FLAGS<m2,s2>);
T_FLAGS<0,0,0> Stp12(V_FLAGS<1>,M_FLAGS<0,0>);
T_FLAGS<0,0,1> Stp12(V_FLAGS<1>,M_FLAGS<1,0>);
template<int m> T_FLAGS<1,1,1> Stp12(V_FLAGS<1>,M_FLAGS<m,1>);

template<class TV> TENSOR<TV> operator+(const TENSOR<TV>& a,const TENSOR<TV>& b)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=a.x(i)+b.x(i);
    return t;
}
template<class TV> TENSOR<TV> operator+(const TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator+(const TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) t.x(i)+=b.v(i);
    return t;
}
template<class TV> TENSOR<TV> operator+(const TENSOR<TV>& a,const VEC_ID_TENSOR_1<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) for(int j=0;j<TV::m;j++) t.x(i)(j,i)+=b.v(j);
    return t;
}
template<class TV> TENSOR<TV> operator+(const TENSOR<TV>& a,const VEC_ID_TENSOR_2<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) for(int j=0;j<TV::m;j++) t.x(i)(i,j)+=b.v(j);
    return t;
}
template<class TV> TENSOR<TV> operator+(const TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) for(int j=0;j<TV::m;j++){t.x(i)(j,i)+=b.v(j);t.x(i)(i,j)+=b.v(j);}
    return t;
}
template<class TV> TENSOR<TV> operator+(const TENSOR<TV>& a,const PERM_TENSOR<TV>& b)
{
    TENSOR<TV> t=a;
    t.x(0)(1,2)+=b.x;
    t.x(1)(2,0)+=b.x;
    t.x(2)(0,1)+=b.x;
    t.x(0)(2,1)-=b.x;
    t.x(1)(0,2)-=b.x;
    t.x(2)(1,0)-=b.x;
    return t;
}

template<class TV> SYMMETRIC_TENSOR<TV> operator+(const SYMMETRIC_TENSOR<TV>& a,const SYMMETRIC_TENSOR<TV>& b)
{
    SYMMETRIC_TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=a.x(i)+b.x(i);
    return t;
}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const SYMMETRIC_TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const SYMMETRIC_TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b)
{
    SYMMETRIC_TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) t.x(i)+=b.v(i);
    return t;
}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const SYMMETRIC_TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b)
{
    SYMMETRIC_TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++){
        for(int j=0;j<=i;j++) t.x(i).Element_Upper(j,i)+=b.v(j);
        for(int j=i;j<TV::m;j++) t.x(i).Element_Lower(j,i)+=b.v(j);}
    return t;
}

template<class TV> TENSOR<TV> operator+(const ZERO_TENSOR<TV>& a,const TENSOR<TV>& b){return b;}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const ZERO_TENSOR<TV>& a,const SYMMETRIC_TENSOR<TV>& b){return b;}
template<class TV> ZERO_TENSOR<TV> operator+(const ZERO_TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> VEC_ID_TENSOR_0<TV> operator+(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return b;}
template<class TV> VEC_ID_TENSOR_1<TV> operator+(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return b;}
template<class TV> VEC_ID_TENSOR_2<TV> operator+(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return b;}
template<class TV> VEC_ID_TENSOR_12<TV> operator+(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return b;}
template<class TV> PERM_TENSOR<TV> operator+(const ZERO_TENSOR<TV>& a,const PERM_TENSOR<TV>& b){return b;}

template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const TENSOR<TV>& b){return b+a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const SYMMETRIC_TENSOR<TV>& b){return b+a;}
template<class TV> VEC_ID_TENSOR_0<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> VEC_ID_TENSOR_0<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return VEC_ID_TENSOR_0<TV>(a.v+b.v);}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return TENSOR<TV>()+a+b;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return TENSOR<TV>()+a+b;}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return SYMMETRIC_TENSOR<TV>()+a+b;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_0<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a+b;}

template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_1<TV>& a,const TENSOR<TV>& b){return b+a;}
template<class TV> VEC_ID_TENSOR_1<TV> operator+(const VEC_ID_TENSOR_1<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return b+a;}
template<class TV> VEC_ID_TENSOR_1<TV> operator+(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return VEC_ID_TENSOR_1<TV>(a.v+b.v);}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return TENSOR<TV>()+a+b;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return TENSOR<TV>()+a+b;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_1<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a+b;}

template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_2<TV>& a,const TENSOR<TV>& b){return b+a;}
template<class TV> VEC_ID_TENSOR_2<TV> operator+(const VEC_ID_TENSOR_2<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return b+a;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return b+a;}
template<class TV> VEC_ID_TENSOR_2<TV> operator+(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return VEC_ID_TENSOR_2<TV>(a.v+b.v);}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return TENSOR<TV>()+a+b;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_2<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a+b;}

template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const TENSOR<TV>& b){return b+a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const SYMMETRIC_TENSOR<TV>& b){return b+a;}
template<class TV> VEC_ID_TENSOR_12<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return b+a;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return b+a;}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return b+a;}
template<class TV> VEC_ID_TENSOR_12<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return VEC_ID_TENSOR_12<TV>(a.v+b.v);}
template<class TV> TENSOR<TV> operator+(const VEC_ID_TENSOR_12<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a+b;}

template<class TV> TENSOR<TV> operator+(const PERM_TENSOR<TV>& a,const TENSOR<TV>& b){return b+a;}
template<class TV> PERM_TENSOR<TV> operator+(const PERM_TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator+(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return b+a;}
template<class TV> TENSOR<TV> operator+(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return b+a;}
template<class TV> TENSOR<TV> operator+(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return b+a;}
template<class TV> TENSOR<TV> operator+(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return b+a;}
template<class TV> PERM_TENSOR<TV> operator+(const PERM_TENSOR<TV>& a,const PERM_TENSOR<TV>& b){return PERM_TENSOR<TV>(a.x+b.x);}

template<class TV> TENSOR<TV> operator-(const TENSOR<TV>& a,const TENSOR<TV>& b)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=a.x(i)-b.x(i);
    return t;
}
template<class TV> TENSOR<TV> operator-(const TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator-(const TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) t.x(i)-=b.v(i);
    return t;
}
template<class TV> TENSOR<TV> operator-(const TENSOR<TV>& a,const VEC_ID_TENSOR_1<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) for(int j=0;j<TV::m;j++) t.x(i)(j,i)-=b.v(j);
    return t;
}
template<class TV> TENSOR<TV> operator-(const TENSOR<TV>& a,const VEC_ID_TENSOR_2<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) for(int j=0;j<TV::m;j++) t.x(i)(i,j)-=b.v(j);
    return t;
}
template<class TV> TENSOR<TV> operator-(const TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b)
{
    TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) for(int j=0;j<TV::m;j++){t.x(i)(j,i)-=b.v(j);t.x(i)(i,j)-=b.v(j);}
    return t;
}
template<class TV> TENSOR<TV> operator-(const TENSOR<TV>& a,const PERM_TENSOR<TV>& b)
{
    TENSOR<TV> t=a;
    t.x(0)(1,2)-=b.x;
    t.x(1)(2,0)-=b.x;
    t.x(2)(0,1)-=b.x;
    t.x(0)(2,1)+=b.x;
    t.x(1)(0,2)+=b.x;
    t.x(2)(1,0)+=b.x;
    return t;
}

template<class TV> SYMMETRIC_TENSOR<TV> operator-(const SYMMETRIC_TENSOR<TV>& a,const SYMMETRIC_TENSOR<TV>& b)
{
    SYMMETRIC_TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=a.x(i)-b.x(i);
    return t;
}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const SYMMETRIC_TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const SYMMETRIC_TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b)
{
    SYMMETRIC_TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++) t.x(i)-=b.v(i);
    return t;
}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const SYMMETRIC_TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b)
{
    SYMMETRIC_TENSOR<TV> t=a;
    for(int i=0;i<TV::m;i++){
        for(int j=0;j<=i;j++) t.x(i).Element_Upper(j,i)-=b.v(j);
        for(int j=i;j<TV::m;j++) t.x(i).Element_Lower(j,i)-=b.v(j);}
    return t;
}

template<class TV> TENSOR<TV> operator-(const ZERO_TENSOR<TV>& a,const TENSOR<TV>& b){return -b;}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const ZERO_TENSOR<TV>& a,const SYMMETRIC_TENSOR<TV>& b){return -b;}
template<class TV> ZERO_TENSOR<TV> operator-(const ZERO_TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> VEC_ID_TENSOR_0<TV> operator-(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return -b;}
template<class TV> VEC_ID_TENSOR_1<TV> operator-(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return -b;}
template<class TV> VEC_ID_TENSOR_2<TV> operator-(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return -b;}
template<class TV> VEC_ID_TENSOR_12<TV> operator-(const ZERO_TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return -b;}
template<class TV> PERM_TENSOR<TV> operator-(const ZERO_TENSOR<TV>& a,const PERM_TENSOR<TV>& b){return -b;}

template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const TENSOR<TV>& b){return -b+a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const SYMMETRIC_TENSOR<TV>& b){return -b+a;}
template<class TV> VEC_ID_TENSOR_0<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> VEC_ID_TENSOR_0<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return VEC_ID_TENSOR_0<TV>(a.v-b.v);}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return SYMMETRIC_TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_0<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a-b;}

template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_1<TV>& a,const TENSOR<TV>& b){return -b+a;}
template<class TV> VEC_ID_TENSOR_1<TV> operator-(const VEC_ID_TENSOR_1<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> VEC_ID_TENSOR_1<TV> operator-(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return VEC_ID_TENSOR_1<TV>(a.v-b.v);}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_1<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_1<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a-b;}

template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_2<TV>& a,const TENSOR<TV>& b){return -b+a;}
template<class TV> VEC_ID_TENSOR_2<TV> operator-(const VEC_ID_TENSOR_2<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> VEC_ID_TENSOR_2<TV> operator-(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return VEC_ID_TENSOR_2<TV>(a.v-b.v);}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_2<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_2<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a-b;}

template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const TENSOR<TV>& b){return -b+a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const SYMMETRIC_TENSOR<TV>& b){return -b+a;}
template<class TV> VEC_ID_TENSOR_12<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> SYMMETRIC_TENSOR<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return SYMMETRIC_TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> VEC_ID_TENSOR_12<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return VEC_ID_TENSOR_12<TV>(a.v-b.v);}
template<class TV> TENSOR<TV> operator-(const VEC_ID_TENSOR_12<TV>& a,const PERM_TENSOR<TV>& b){return TENSOR<TV>()+a-b;}

template<class TV> TENSOR<TV> operator-(const PERM_TENSOR<TV>& a,const TENSOR<TV>& b){return -b+a;}
template<class TV> PERM_TENSOR<TV> operator-(const PERM_TENSOR<TV>& a,const ZERO_TENSOR<TV>& b){return a;}
template<class TV> TENSOR<TV> operator-(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_0<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_1<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_2<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> TENSOR<TV> operator-(const PERM_TENSOR<TV>& a,const VEC_ID_TENSOR_12<TV>& b){return TENSOR<TV>()+a-b;}
template<class TV> PERM_TENSOR<TV> operator-(const PERM_TENSOR<TV>& a,const PERM_TENSOR<TV>& b){return PERM_TENSOR<TV>(a.x-b.x);}

template<class TV> SYMMETRIC_TENSOR<TV> Contract_0(const SYMMETRIC_TENSOR<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    SYMMETRIC_TENSOR<TV> t;
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            t.x(j)+=a.x(i)*m(i,j);
    return t;
}
template<class TV> TENSOR<TV> Contract_0(const TENSOR<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            t.x(j)+=a.x(i)*m(i,j);
    return t;
}
template<class TV> ZERO_TENSOR<TV> Contract_0(const ZERO_TENSOR<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m){return a;}
template<class TV> VEC_ID_TENSOR_0<TV> Contract_0(const VEC_ID_TENSOR_0<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m){return VEC_ID_TENSOR_0<TV>(m.Transpose_Times(a.v));}
template<class TV> TENSOR<TV> Contract_0(const VEC_ID_TENSOR_1<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m){return Tensor_Product_1(m.Transposed(),a.v);}
template<class TV> TENSOR<TV> Contract_0(const VEC_ID_TENSOR_2<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m){return Tensor_Product_2(m.Transposed(),a.v);}
template<class TV> SYMMETRIC_TENSOR<TV> Contract_0(const VEC_ID_TENSOR_12<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m){return Symmetric_Tensor_Product_12(m.Transposed(),a.v);}
template<class TV> TENSOR<TV> Contract_0(const PERM_TENSOR<TV>& a,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=MATRIX<typename TV::SCALAR,TV::m>::Cross_Product_Matrix(-a.x*m.Column(i));
    return t;
}

template<int m,int a,int b> T_FLAGS<0,0,0> C0(T_FLAGS<m,a,b>,M_FLAGS<0,0>);
template<int m,int a,int b> T_FLAGS<m,a,b> C0(T_FLAGS<m,a,b>,M_FLAGS<0,1>);
template<int m,int a,int b> T_FLAGS<m,a,b> C0(T_FLAGS<m,a,b>,M_FLAGS<1,1>);
T_FLAGS<0,0,0> C0(T_FLAGS<0,0,0>,M_FLAGS<1,0>);
T_FLAGS<1,0,0> C0(T_FLAGS<1,0,0>,M_FLAGS<1,0>);
template<int m,int a,int b> T_FLAGS<0,0,1> C0(T_FLAGS<m,a,b>,M_FLAGS<1,0>);

template<class TV> ZERO_TENSOR<TV> Contract_1(const ZERO_TENSOR<TV>& p,const MATRIX<typename TV::SCALAR,TV::m>& m) {return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Contract_1(const ZERO_TENSOR<TV>& p,const ZERO_MAT<TV>& m) {return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Contract_1(const ZERO_TENSOR<TV>& p,ID_MAT<TV> m) {return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Contract_1(const ZERO_TENSOR<TV>& p,SCALE_MAT<TV> m) {return ZERO_TENSOR<TV>();}

template<class TV> TENSOR<TV> Contract_1(const TENSOR<TV>& p,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=m.Transpose_Times(p.x(i));
    return t;
}
template<class TV> ZERO_TENSOR<TV> Contract_1(const TENSOR<TV>& p,const ZERO_MAT<TV>& m) {return ZERO_TENSOR<TV>();}
template<class TV> TENSOR<TV> Contract_1(const TENSOR<TV>& p,ID_MAT<TV> m) {return p;}
template<class TV> TENSOR<TV> Contract_1(const TENSOR<TV>& p,SCALE_MAT<TV> m) {return p*m.x;}

template<class TV> TENSOR<TV> Contract_1(const PERM_TENSOR<TV>& p,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++){
        MATRIX<typename TV::SCALAR,TV::m> c=MATRIX<typename TV::SCALAR,TV::m>::Cross_Product_Matrix(p.x*m.Column(i));
        for(int j=0;j<TV::m;j++) t.x(j).Set_Row(i,c.Row(j));}
    return t;
}
template<class TV> ZERO_TENSOR<TV> Contract_1(const PERM_TENSOR<TV>& p,const ZERO_MAT<TV>& m) {return ZERO_TENSOR<TV>();}
template<class TV> PERM_TENSOR<TV> Contract_1(const PERM_TENSOR<TV>& p,ID_MAT<TV> m) {return p;}
template<class TV> PERM_TENSOR<TV> Contract_1(const PERM_TENSOR<TV>& p,SCALE_MAT<TV> m) {return p*m.x;}

template<int m,int a,int b> T_FLAGS<0,0,0> C1(T_FLAGS<m,a,b>,M_FLAGS<0,0>);
template<int m,int a,int b> T_FLAGS<m,a,b> C1(T_FLAGS<m,a,b>,M_FLAGS<0,1>);
template<int m,int a,int b> T_FLAGS<m,a,b> C1(T_FLAGS<m,a,b>,M_FLAGS<1,1>);
T_FLAGS<0,0,0> C1(T_FLAGS<0,0,0>,M_FLAGS<1,0>);
T_FLAGS<1,0,1> C1(T_FLAGS<1,0,1>,M_FLAGS<1,0>);
template<int m,int a,int b> T_FLAGS<0,0,1> C1(T_FLAGS<m,a,b>,M_FLAGS<1,0>);

template<class TV> TENSOR<TV> Contract_2(const TENSOR<TV>& p,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++) t.x(i)=p.x(i)*m;
    return t;
}
template<class TV> TENSOR<TV> Contract_2(const PERM_TENSOR<TV>& p,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    TENSOR<TV> t;
    for(int i=0;i<TV::m;i++){
        MATRIX<typename TV::SCALAR,TV::m> c=MATRIX<typename TV::SCALAR,TV::m>::Cross_Product_Matrix(p.x*m.Column(i));
        for(int j=0;j<TV::m;j++) t.x(j).Set_Column(i,c.Column(j));}
    return t;
}
template<class TV> ZERO_TENSOR<TV> Contract_2(const PERM_TENSOR<TV>& p,const ZERO_MAT<TV>& m) {return ZERO_TENSOR<TV>();}
template<class TV> PERM_TENSOR<TV> Contract_2(const PERM_TENSOR<TV>& p,ID_MAT<TV> m) {return p;}
template<class TV> PERM_TENSOR<TV> Contract_2(const PERM_TENSOR<TV>& p,SCALE_MAT<TV> m) {return p*m.x;}

template<int m,int a,int b> T_FLAGS<0,0,0> C2(T_FLAGS<m,a,b>,M_FLAGS<0,0>);
template<int m,int a,int b> T_FLAGS<m,a,b> C2(T_FLAGS<m,a,b>,M_FLAGS<0,1>);
template<int m,int a,int b> T_FLAGS<m,a,b> C2(T_FLAGS<m,a,b>,M_FLAGS<1,1>);
T_FLAGS<0,0,0> C2(T_FLAGS<0,0,0>,M_FLAGS<1,0>);
T_FLAGS<1,1,0> C2(T_FLAGS<1,1,0>,M_FLAGS<1,0>);
template<int m,int a,int b> T_FLAGS<0,0,1> C2(T_FLAGS<m,a,b>,M_FLAGS<1,0>);

template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ZERO_MAT<TV>& cm1,const ZERO_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ZERO_MAT<TV>& cm1,const ID_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ZERO_MAT<TV>& cm1,const SCALE_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ZERO_MAT<TV>& cm1,const MATRIX<typename TV::SCALAR,TV::m>){return ZERO_TENSOR<TV>();}

template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ID_MAT<TV>& cm1,const ZERO_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ID_MAT<TV>& cm1,const ID_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ID_MAT<TV>& cm1,const SCALE_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> SYMMETRIC_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const ID_MAT<TV>& cm1,const MATRIX<typename TV::SCALAR,TV::m>& cm2)
{return Contract_2(t,cm2).Twice_Symmetric_Part_12();}

template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const SCALE_MAT<TV>& cm1,const ZERO_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const SCALE_MAT<TV>& cm1,const ID_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const SCALE_MAT<TV>& cm1,const SCALE_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> SYMMETRIC_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const SCALE_MAT<TV>& cm1,const MATRIX<typename TV::SCALAR,TV::m>& cm2)
{return Symmetric_Double_Contract_12(t*cm1.x,ID_MAT<TV>(),cm2);}

template<class TV> ZERO_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const MATRIX<typename TV::SCALAR,TV::m>& cm1,const ZERO_MAT<TV>& cm2){return ZERO_TENSOR<TV>();}
template<class TV> SYMMETRIC_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const MATRIX<typename TV::SCALAR,TV::m>& cm1,const ID_MAT<TV>& cm2)
{return Symmetric_Double_Contract_12(-t,cm2,cm1);}
template<class TV> SYMMETRIC_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const MATRIX<typename TV::SCALAR,TV::m>& cm1,const SCALE_MAT<TV>& cm2)
{return Symmetric_Double_Contract_12(-t,cm2,cm1);}
template<class TV> SYMMETRIC_TENSOR<TV> Symmetric_Double_Contract_12(const PERM_TENSOR<TV>&t,const MATRIX<typename TV::SCALAR,TV::m>& cm1,const MATRIX<typename TV::SCALAR,TV::m>& cm2)
{
    SYMMETRIC_TENSOR<TV> s;
    TV c10=cm1.Row(0),c11=cm1.Row(1),c12=cm1.Row(2),c20=t.x*cm2.Row(0),c21=t.x*cm2.Row(1),c22=t.x*cm2.Row(2);
    s.x(0)=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(c11,c22)-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(c12,c21);
    s.x(1)=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(c12,c20)-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(c10,c22);
    s.x(2)=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(c10,c21)-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(c11,c20);
    return s;
}

template<int m,int s,int m2,int s2> T_FLAGS<0,0,(m&~s&(m2|s2))|(m2&~s2&(m|s))> Sdc12(T_FLAGS<0,1,1>,M_FLAGS<m,s>,M_FLAGS<m2,s2>);
}
}
#endif
