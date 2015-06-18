//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_DIFF_EXT
//##################################################################### 
#ifndef __AUTO_DIFF_EXT__
#define __AUTO_DIFF_EXT__

#include <Tools/Auto_Diff/GRADIENT.h>
#include <Tools/Auto_Diff/GRADIENT_VEC.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
using ::std::sin;
using ::std::cos;
using ::std::tan;
using ::std::exp;
using ::std::log;
using ::std::abs;
using ::std::sqrt;
using ::PhysBAM::sqr;
using ::PhysBAM::cube;

template<class T,class VEC> struct AUTO_DIFF_EXT;

template<class T,class VEC> AUTO_DIFF_EXT<T,VEC>
Make_Diff(T x,const GRADIENT<T,VEC>& dx);

template<class T,class VEC>
struct AUTO_DIFF_EXT
{
    T x;
    GRADIENT<T,VEC> dx;

    AUTO_DIFF_EXT(T z=T()): x(z)
    {}

    AUTO_DIFF_EXT(T z,const GRADIENT<T,VEC>& dz)
        :x(z),dx(dz)
    {}

    decltype(Make_Diff(-x,-dx)) operator-() const
    {return Make_Diff(-x,-dx);}

    AUTO_DIFF_EXT operator+() const
    {return *this;}

    template<class VEC1> auto
    operator+(const AUTO_DIFF_EXT<T,VEC1>& a) const -> decltype(Make_Diff(x+a.x,dx+a.dx))
    {return Make_Diff(x+a.x,dx+a.dx);}

    template<class VEC1> auto
    operator-(const AUTO_DIFF_EXT<T,VEC1>& a) const -> decltype(Make_Diff(x-a.x,dx-a.dx))
    {return Make_Diff(x-a.x,dx-a.dx);}

    template<class VEC1> auto
    operator*(const AUTO_DIFF_EXT<T,VEC1>& a) const -> decltype(Make_Diff(x*a.x,a.x*dx+x*a.dx))
    {return Make_Diff(x*a.x,a.x*dx+x*a.dx);}

    template<class VEC1> auto
    operator/(const AUTO_DIFF_EXT<T,VEC1>& a) const -> decltype(Make_Diff(T(),(dx-T()*a.dx)/a.x))
    {T z=x/a.x;auto q=dx-z*a.dx;return Make_Diff(z,q/a.x);}

    AUTO_DIFF_EXT operator+(T a) const
    {return Make_Diff(x+a,dx);}

    AUTO_DIFF_EXT operator-(T a) const
    {return Make_Diff(x-a,dx);}

    auto operator*(T a) const -> decltype(Make_Diff(x*a,a*dx))
    {return Make_Diff(x*a,a*dx);}

    auto operator/(T a) const -> decltype(Make_Diff(x/a,dx/a))
    {return Make_Diff(x/a,dx/a);}

    template<class VEC1>
    void Fill_From(const AUTO_DIFF_EXT<T,VEC1>& z)
    {x=z.x;::PhysBAM::HETERO_DIFF::Fill_From(dx,z.dx);}
};

template<class T,class VEC> AUTO_DIFF_EXT<T,VEC>
Make_Diff(T x,const GRADIENT<T,VEC>& dx)
{return AUTO_DIFF_EXT<T,VEC>(x,dx);}

template<class T,class LAYOUT>
AUTO_DIFF_EXT<T,typename EMPTY_VEC<LAYOUT,-1>::TYPE>
Diff_From_Const(T a)
{return AUTO_DIFF_EXT<T,typename EMPTY_VEC<LAYOUT,-1>::TYPE>(a);}

template<class T,class VEC>
inline auto operator+(T c,const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(c+a.x,a.dx))
{return Make_Diff(c+a.x,a.dx);}

template<class T,class VEC>
inline auto operator-(T c,const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(c-a.x,-a.dx))
{return Make_Diff(c-a.x,-a.dx);}

template<class T,class VEC> inline auto
operator*(T c,const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(c*a.x,c*a.dx))
{return Make_Diff(c*a.x,c*a.dx);}

template<class T,class VEC> inline auto
operator/(T c,const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(T(),-T()*a.dx))
{T z=c/a.x,w=z/a.x;return Make_Diff(z,-w*a.dx);}

template<class T,class VEC> inline auto
sqrt(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(T(),a.dx/(2*T())))
{T s=sqrt(a.x);auto t=a.dx/(2*s);return Make_Diff(s,t);}

template<class T,class VEC> inline auto
sqr(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(sqr(a.x),2*a.x*a.dx))
{return Make_Diff(sqr(a.x),2*a.x*a.dx);}

template<class T,class VEC> inline auto
cube(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(a.x*T(),3*T()*a.dx))
{T sq=sqr(a.x);return Make_Diff(a.x*sq,3*sq*a.dx);}

template<class T,class VEC,class VEC1> inline AUTO_DIFF_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))>
max(const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC1>& b)
{
    AUTO_DIFF_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))> r;
    if(a.x>b.x) r.Fill_From(a);
    else r.Fill_From(b);
    return r;
}

template<class T,class VEC> inline decltype(AUTO_DIFF_EXT<T,VEC>()+T())
max(const AUTO_DIFF_EXT<T,VEC>& a,T b)
{
    decltype(AUTO_DIFF_EXT<T,VEC>()+T()) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class T,class VEC> inline decltype(AUTO_DIFF_EXT<T,VEC>()+T())
max(T b,const AUTO_DIFF_EXT<T,VEC>& a)
{
    decltype(AUTO_DIFF_EXT<T,VEC>()+T()) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class T,class VEC>
inline AUTO_DIFF_EXT<T,VEC>
max(const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC>& b)
{return a.x>b.x?a:b;}

template<class T,class VEC,class VEC1> inline AUTO_DIFF_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))>
min(const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC1>& b)
{
    AUTO_DIFF_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))> r;
    if(a.x>b.x) r.Fill_From(b);
    else r.Fill_From(a);
    return r;
}

template<class T,class VEC> inline decltype(AUTO_DIFF_EXT<T,VEC>()+T())
min(const AUTO_DIFF_EXT<T,VEC>& a,T b)
{
    decltype(AUTO_DIFF_EXT<T,VEC>()+T()) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class T,class VEC> inline decltype(AUTO_DIFF_EXT<T,VEC>()+T())
min(T b,const AUTO_DIFF_EXT<T,VEC>& a)
{
    decltype(AUTO_DIFF_EXT<T,VEC>()+T()) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class T,class VEC>
inline AUTO_DIFF_EXT<T,VEC>
min(const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC>& b)
{return a.x>b.x?b:a;}

template<class T,class VEC> inline auto
log(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(::std::log(a.x),a.dx/a.x))
{auto z=a.dx/a.x;return Make_Diff(::std::log(a.x),z);}

template<class T,class VEC> inline auto
exp(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(T(),T()*a.dx))
{T s=exp(a.x);return Make_Diff(s,s*a.dx);}

template<class T,class VEC> inline auto
sin(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(T(),T()*a.dx))
{T s=sin(a.x),c=cos(a.x);return Make_Diff(s,c*a.dx);}

template<class T,class VEC> inline auto
cos(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(T(),-T()*a.dx))
{T s=sin(a.x),c=cos(a.x);return Make_Diff(c,-s*a.dx);}

template<class T,class VEC> inline auto
tan(const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(Make_Diff(T(),T()*a.dx))
{T t=tan(a.x),s=1+t*t;return Make_Diff(t,s*a.dx);}

template<class T,class VEC> inline AUTO_DIFF_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC())))>
abs(const AUTO_DIFF_EXT<T,VEC>& a)
{
    AUTO_DIFF_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC())))> r;
    if(a.x>=0) r.Fill_From(a);
    else r.Fill_From(-a);
    return r;
}

template<class T,class VEC,class VEC1> inline auto
hypot(const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC1>& b)
    -> decltype(Make_Diff(T(),T()*a.dx+T()*b.dx))
{
    T c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return Make_Diff(c,d*a.dx+e*b.dx);
}

template<class T,class VEC> inline auto
hypot(const AUTO_DIFF_EXT<T,VEC>& a,T b)
    -> decltype(Make_Diff(T(),T()*a.dx))
{
    T c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c;
    return Make_Diff(c,d*a.dx);
}

template<class T,class VEC> inline auto
hypot(T b,const AUTO_DIFF_EXT<T,VEC>& a)
    -> decltype(hypot(a,b))
{return hypot(a,b);}

template<class T,class VEC,class VEC1,class VEC2> inline auto
hypot(const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC1>& b,const AUTO_DIFF_EXT<T,VEC2>& c)
    -> decltype(Make_Diff(T(),T()*a.dx+T()*b.dx+T()*T().dx))
{
    T s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    return Make_Diff(s,aa*a.dx+bb*b.dx+cc*c.dx);
}

template<class T,class VEC,class VEC1> inline auto
hypot(const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC1>& b,T c)
    -> decltype(Make_Diff(T(),T()*a.dx+T()*b.dx))
{
    T s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s;
    return Make_Diff(s,aa*a.dx+bb*b.dx);
}

template<class T,class VEC,class VEC1> inline auto
hypot(const AUTO_DIFF_EXT<T,VEC>& a,T c,const AUTO_DIFF_EXT<T,VEC1>& b)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class VEC1> inline auto
hypot(T c,const AUTO_DIFF_EXT<T,VEC>& a,const AUTO_DIFF_EXT<T,VEC1>& b)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC> inline auto
hypot(const AUTO_DIFF_EXT<T,VEC>& a,T b,T c)
    -> decltype(Make_Diff(T(),T()*a.dx))
{
    T t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s;
    return Make_Diff(s,d*a.dx);
}

template<class T,class VEC> inline auto
hypot(T b,const AUTO_DIFF_EXT<T,VEC>& a,T c)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC> inline auto
hypot(T b,T c,const AUTO_DIFF_EXT<T,VEC>& a)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class VEC1> inline auto
atan2(const AUTO_DIFF_EXT<T,VEC>& y,const AUTO_DIFF_EXT<T,VEC1>& x)
    -> decltype(Make_Diff(::std::atan2(y.x,x.x),T()*y.dx-T()*x.dx))
{
    T c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    auto f=d*y.dx-e*x.dx;
    return Make_Diff(::std::atan2(y.x,x.x),f);
}

template<class T,class VEC> inline auto
atan2(const AUTO_DIFF_EXT<T,VEC>& y,T x)
    -> decltype(Make_Diff(::std::atan2(y.x,x),T()*y.dx))
{
    T c=sqr(x)+sqr(y.x),d=x/c;
    auto f=d*y.dx;
    return Make_Diff(::std::atan2(y.x,x),f);
}

template<class T,class VEC> inline auto
atan2(T y,const AUTO_DIFF_EXT<T,VEC>& x)
    -> decltype(Make_Diff(::std::atan2(y,x.x),-T()*x.dx))
{
    T c=sqr(x.x)+sqr(y),e=y/c;
    auto f=-e*x.dx;
    return Make_Diff(::std::atan2(y,x.x),f);
}
template<class T,class VEC> struct AUTO_DIFF_EXT_VEC;

template<class TV,class VEC> AUTO_DIFF_EXT_VEC<TV,VEC>
Make_Diff_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx);

template<class TV,class VEC> struct AUTO_DIFF_EXT_VEC;

template<class T,class TV,class VEC1,class VEC2> auto
Cross_Helper(const AUTO_DIFF_EXT_VEC<VECTOR<T,3>,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC2>& b)
    -> decltype(Make_Diff_Vec(a.x.Cross(b.x),MATRIX<T,TV::m>()*b.dx-MATRIX<T,TV::m>()*a.dx));

template<class TV,class VEC1,class TYPE>
typename enable_if<TV::m!=3,const AUTO_DIFF_EXT_VEC<TV,VEC1> >::type
Cross_Helper(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a,TYPE);

template<class TV,class VEC>
struct AUTO_DIFF_EXT_VEC
{
    typedef typename TV::SCALAR T;

    TV x;
    GRADIENT_VEC<TV,VEC> dx;

    AUTO_DIFF_EXT_VEC(TV x=TV()):
        x(x)
    {}

    AUTO_DIFF_EXT_VEC(TV x,const GRADIENT_VEC<TV,VEC>& dx):
        x(x),dx(dx)
    {}

    decltype(Make_Diff_Vec(-TV(),-GRADIENT_VEC<TV,VEC>())) operator-() const
    {return Make_Diff_Vec(-x,-dx);}

    AUTO_DIFF_EXT_VEC operator+() const
    {return *this;}

    template<class VEC1>
    auto operator+(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const -> decltype(Make_Diff_Vec(x+a.x,dx+a.dx))
    {return Make_Diff_Vec(x+a.x,dx+a.dx);}

    template<class VEC1>
    auto operator-(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const -> decltype(Make_Diff_Vec(x-a.x,dx-a.dx))
    {return Make_Diff_Vec(x-a.x,dx-a.dx);}

    template<class VEC1> auto
    operator*(const AUTO_DIFF_EXT<T,VEC1>& a) const -> decltype(Make_Diff_Vec(x*a.x,a.x*dx+Outer_Product(x,a.dx)))
    {return Make_Diff_Vec(x*a.x,a.x*dx+Outer_Product(x,a.dx));}

    template<class VEC1> auto
    operator/(const AUTO_DIFF_EXT<T,VEC1>& a) const -> decltype(Make_Diff_Vec(x/a.x,dx/a.x+Outer_Product(x,-a.dx/(a.x*a.x))))
    {auto p=-a.dx/(a.x*a.x);return Make_Diff_Vec(x/a.x,dx/a.x+Outer_Product(x,p));}

    AUTO_DIFF_EXT_VEC operator+(TV a) const
    {return Make_Diff_Vec(x+a,dx);}

    AUTO_DIFF_EXT_VEC operator-(TV a) const
    {return Make_Diff_Vec(x-a,dx);}

    auto operator*(T a) const -> decltype(Make_Diff_Vec(x*a,dx*a))
    {return Make_Diff_Vec(x*a,dx*a);}

    auto operator/(T a) const -> decltype(Make_Diff_Vec(x/a,dx/a))
    {return Make_Diff_Vec(x/a,dx/a);}

    auto Dot(TV v) const -> decltype(Make_Diff(x.Dot(v),dx.Transpose_Times(v)))
    {return Make_Diff(x.Dot(v),dx.Transpose_Times(v));}

    template<class VEC1> auto
    Dot(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const
        -> decltype(Make_Diff(x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x)))
    {return Make_Diff(x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x));}

    auto Cross(const TV& a) const -> decltype(Make_Diff_Vec(x.Cross(a),-MATRIX<T,TV::m>()*dx))
    {
        MATRIX<T,TV::m> cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(a);
        return Make_Diff_Vec(x.Cross(a),-cp_a*dx);
    }

    template<class VEC1>
    auto Cross(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const -> decltype(Cross_Helper(*this,AUTO_DIFF_EXT_VEC<TV,VEC1>()))
    {return Cross_Helper(*this,a);}

    auto Magnitude_Squared() const -> decltype(Make_Diff(x.Magnitude_Squared(),dx.Transpose_Times(x)*2))
    {return Make_Diff(x.Magnitude_Squared(),dx.Transpose_Times(x)*2);}

    auto Magnitude() const -> decltype(Make_Diff(T(),dx.Transpose_Times(x)/T()))
    {T s=x.Magnitude();auto t=dx.Transpose_Times(x)/s;return Make_Diff(s,t);}
};


template<class T,class TV,class VEC1,class VEC2> auto
Cross_Helper(const AUTO_DIFF_EXT_VEC<VECTOR<T,3>,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC2>& b)
    -> decltype(Make_Diff_Vec(a.x.Cross(b.x),MATRIX<T,TV::m>()*b.dx-MATRIX<T,TV::m>()*a.dx))
{
    MATRIX<T,TV::m> cp_t=MATRIX<T,TV::m>::Cross_Product_Matrix(a.x),cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(b.x);
    return Make_Diff_Vec(a.x.Cross(b.x),cp_t*b.dx-cp_a*a.dx);
}

template<class TV,class VEC1,class TYPE>
typename enable_if<TV::m!=3,const AUTO_DIFF_EXT_VEC<TV,VEC1> >::type
Cross_Helper(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a,TYPE)
{
    PHYSBAM_FATAL_ERROR("Cross product not defined except in 3D.");
    return a;
}

template<class TV,class VEC> inline AUTO_DIFF_EXT_VEC<TV,VEC>
Make_Diff_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx)
{return AUTO_DIFF_EXT_VEC<TV,VEC>(x,dx);}

template<class LAYOUT,int i,class TV> inline
AUTO_DIFF_EXT_VEC<TV,typename ONE_NONZERO_VECTOR<i,LAYOUT,TV::m>::TYPE> Diff_From_Var(const TV& v)
{return AUTO_DIFF_EXT_VEC<TV,typename ONE_NONZERO_VECTOR<i,LAYOUT,TV::m>::TYPE>(v);}

template<class T,int d,class VEC> inline auto
operator*(const AUTO_DIFF_EXT<T,VEC>& a,const VECTOR<T,d>& v) -> decltype(Make_Diff_Vec(v*a.x,Outer_Product(v,a.dx)))
{return Make_Diff_Vec(v*a.x,Outer_Product(v,a.dx));}

template<class T,int d,class VEC> inline auto
operator*(const VECTOR<T,d>& v,const AUTO_DIFF_EXT<T,VEC>& a) -> decltype(a*v)
{return a*v;}

template<class T,class TV,class VEC,class VEC1> inline auto
operator*(const AUTO_DIFF_EXT<T,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC>& v) -> decltype(v*a)
{return v*a;}

template<class TV,class VEC> auto
operator*(typename TV::SCALAR a,const AUTO_DIFF_EXT_VEC<TV,VEC>& v) -> decltype(v*a)
{return v*a;}

template<class TV,class VEC> auto
operator+(TV a,const AUTO_DIFF_EXT_VEC<TV,VEC>& u) -> decltype(Make_Diff_Vec(u.x+a,u.dx))
{return Make_Diff_Vec(u.x+a,u.dx);}

template<class TV,class VEC> auto
operator-(TV a,const AUTO_DIFF_EXT_VEC<TV,VEC>& u) -> decltype(Make_Diff_Vec(a-u.x,-u.dx))
{return Make_Diff_Vec(a-u.x,-u.dx);}

template<class T,int d,class VEC1> auto
operator/(const VECTOR<T,d>& u,const AUTO_DIFF_EXT<T,VEC1>& a) -> decltype(u*((T)1/a))
{return u*((T)1/a);}

template<class TV,class VEC> auto
operator*(const MATRIX<typename TV::SCALAR,TV::m>& m,const AUTO_DIFF_EXT_VEC<TV,VEC>& v)
    -> decltype(Make_Diff_Vec(m*v.x,m*v.dx))
{return Make_Diff_Vec(m*v.x,m*v.dx);}
}
using HETERO_DIFF::AUTO_DIFF_EXT;
using HETERO_DIFF::AUTO_DIFF_EXT_VEC;
using HETERO_DIFF::Diff_From_Const;
using HETERO_DIFF::Diff_From_Var;
using HETERO_DIFF::Extract;
using HETERO_DIFF::sin;
using HETERO_DIFF::cos;
using HETERO_DIFF::tan;
//using HETERO_DIFF::atan2;
using HETERO_DIFF::exp;
using HETERO_DIFF::log;
// using HETERO_DIFF::min;
// using HETERO_DIFF::max;
using HETERO_DIFF::abs;
//using HETERO_DIFF::hypot;
using HETERO_DIFF::sqr;
using HETERO_DIFF::sqrt;
}
#endif
