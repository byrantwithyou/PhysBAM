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

#define SC typename TV::SCALAR()
#define DX GRADIENT<TV,VEC>()
#define DX1 GRADIENT<TV,VEC1>()
#define DX2 GRADIENT<TV,VEC2>()
#define VDX GRADIENT_VEC<TV,VEC>()
#define VDX1 GRADIENT_VEC<TV,VEC1>()
#define VDX2 GRADIENT_VEC<TV,VEC2>()
template<class TV,class VEC> struct AUTO_DIFF_EXT;

template<class TV,class VEC> AUTO_DIFF_EXT<TV,VEC>
Make_Diff(typename TV::SCALAR x,const GRADIENT<TV,VEC>& dx);

template<class TV,class VEC>
struct AUTO_DIFF_EXT
{
    typedef typename TV::SCALAR T;

    T x;
    GRADIENT<TV,VEC> dx;

    AUTO_DIFF_EXT(T z=T()): x(z)
    {}

    AUTO_DIFF_EXT(T z,const GRADIENT<TV,VEC>& dz)
        :x(z),dx(dz)
    {}

    decltype(Make_Diff(-x,-dx)) operator-() const
    {return Make_Diff(-x,-dx);}

    AUTO_DIFF_EXT operator+() const
    {return *this;}

    template<class VEC1> decltype(Make_Diff(x+SC,dx+DX1))
    operator+(const AUTO_DIFF_EXT<TV,VEC1>& a) const
    {return Make_Diff(x+a.x,dx+a.dx);}

    template<class VEC1> decltype(Make_Diff(x-SC,dx-DX1))
    operator-(const AUTO_DIFF_EXT<TV,VEC1>& a) const
    {return Make_Diff(x-a.x,dx-a.dx);}

    template<class VEC1> decltype(Make_Diff(x*SC,SC*dx+x*DX1))
    operator*(const AUTO_DIFF_EXT<TV,VEC1>& a) const
    {return Make_Diff(x*a.x,a.x*dx+x*a.dx);}

    template<class VEC1> decltype(Make_Diff(x/SC,(dx-(x/SC)*DX1)/SC))
    operator/(const AUTO_DIFF_EXT<TV,VEC1>& a) const
    {T z=x/a.x,ax2=sqr(a.x);auto q=dx-z*a.dx;return Make_Diff(z,q/a.x);}

    AUTO_DIFF_EXT operator+(T a) const
    {return Make_Diff(x+a,dx);}

    AUTO_DIFF_EXT operator-(T a) const
    {return Make_Diff(x-a,dx);}

    decltype(Make_Diff(x*SC,dx*SC)) operator*(T a) const
    {return Make_Diff(x*a,a*dx);}

    decltype(Make_Diff(x/SC,dx/SC)) operator/(T a) const
    {return Make_Diff(x/a,dx/a);}

    template<class VEC1>
    void Fill_From(const AUTO_DIFF_EXT<TV,VEC1>& z)
    {x=z.x;::PhysBAM::HETERO_DIFF::Fill_From(dx,z.dx);}
};

template<class TV,class VEC> AUTO_DIFF_EXT<TV,VEC>
Make_Diff(typename TV::SCALAR x,const GRADIENT<TV,VEC>& dx)
{return AUTO_DIFF_EXT<TV,VEC>(x,dx);}

template<int n,int i,class TV>
AUTO_DIFF_EXT<TV,typename ONE_NONZERO_VECTOR<TV,n,i>::TYPE>
Diff_From_Var(TV v,int j)
{
    AUTO_DIFF_EXT<TV,typename ONE_NONZERO_VECTOR<TV,n,i>::TYPE> r(v(j));
    Set<i>(r.dx.x,TV::Axis_Vector(j));
    return r;
}

template<class TV,int n>
AUTO_DIFF_EXT<TV,typename EMPTY_VEC<TV,n>::TYPE>
Diff_From_Const(typename TV::SCALAR a)
{return AUTO_DIFF_EXT<TV,typename EMPTY_VEC<TV,n>::TYPE>(a);}

template<class TV,class VEC>
inline AUTO_DIFF_EXT<TV,VEC> operator+(typename TV::SCALAR c,const AUTO_DIFF_EXT<TV,VEC>& a)
{return Make_Diff(c+a.x,a.dx);}

template<class TV,class VEC>
inline AUTO_DIFF_EXT<TV,VEC> operator-(typename TV::SCALAR c,const AUTO_DIFF_EXT<TV,VEC>& a)
{return Make_Diff(c-a.x,-a.dx);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
operator*(typename TV::SCALAR c,const AUTO_DIFF_EXT<TV,VEC>& a)
{return Make_Diff(c*a.x,c*a.dx);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
operator/(typename TV::SCALAR c,const AUTO_DIFF_EXT<TV,VEC>& a)
{typename TV::SCALAR z=c/a.x,w=z/a.x;return Make_Diff(z,-w*a.dx);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,DX/SC))
sqrt(const AUTO_DIFF_EXT<TV,VEC>& a)
{typename TV::SCALAR s=sqrt(a.x);auto t=a.dx/(2*s);return Make_Diff(s,t);}

template<class TV,class VEC> inline decltype(Make_Diff(sqr(SC),2*SC*DX))
sqr(const AUTO_DIFF_EXT<TV,VEC>& a)
{return Make_Diff(sqr(a.x),2*a.x*a.dx);}

template<class TV,class VEC> inline decltype(Make_Diff(cube(SC),SC*DX))
cube(const AUTO_DIFF_EXT<TV,VEC>& a)
{typename TV::SCALAR sq=sqr(a.x);return Make_Diff(a.x*sq,3*sq*a.dx);}

template<class TV,class VEC,class VEC1> inline AUTO_DIFF_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))>
max(const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC1>& b)
{
    AUTO_DIFF_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))> r;
    if(a.x>b.x) r.Fill_From(a);
    else r.Fill_From(b);
    return r;
}

template<class TV,class VEC> inline decltype(AUTO_DIFF_EXT<TV,VEC>()*SC)
max(const AUTO_DIFF_EXT<TV,VEC>& a,typename TV::SCALAR b)
{
    decltype(AUTO_DIFF_EXT<TV,VEC>()*SC) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class TV,class VEC> inline decltype(AUTO_DIFF_EXT<TV,VEC>()*SC)
max(typename TV::SCALAR b,const AUTO_DIFF_EXT<TV,VEC>& a)
{
    decltype(AUTO_DIFF_EXT<TV,VEC>()*SC) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class TV,class VEC>
inline AUTO_DIFF_EXT<TV,VEC>
max(const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC>& b)
{return a.x>b.x?a:b;}

template<class TV,class VEC,class VEC1> inline AUTO_DIFF_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))>
min(const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC1>& b)
{
    AUTO_DIFF_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1()))> r;
    if(a.x>b.x) r.Fill_From(b);
    else r.Fill_From(a);
    return r;
}

template<class TV,class VEC> inline decltype(AUTO_DIFF_EXT<TV,VEC>()*SC)
min(const AUTO_DIFF_EXT<TV,VEC>& a,typename TV::SCALAR b)
{
    decltype(AUTO_DIFF_EXT<TV,VEC>()*SC) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class TV,class VEC> inline decltype(AUTO_DIFF_EXT<TV,VEC>()*SC)
min(typename TV::SCALAR b,const AUTO_DIFF_EXT<TV,VEC>& a)
{
    decltype(AUTO_DIFF_EXT<TV,VEC>()*SC) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class TV,class VEC>
inline AUTO_DIFF_EXT<TV,VEC>
min(const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC>& b)
{return a.x>b.x?b:a;}

template<class TV,class VEC> inline decltype(Make_Diff(SC,DX/SC))
log(const AUTO_DIFF_EXT<TV,VEC>& a)
{auto z=a.dx/a.x;return Make_Diff(::std::log(a.x),z);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
exp(const AUTO_DIFF_EXT<TV,VEC>& a)
{typename TV::SCALAR s=exp(a.x);return Make_Diff(s,s*a.dx);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
sin(const AUTO_DIFF_EXT<TV,VEC>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return Make_Diff(s,c*a.dx);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,-SC*DX))
cos(const AUTO_DIFF_EXT<TV,VEC>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return Make_Diff(c,-s*a.dx);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
tan(const AUTO_DIFF_EXT<TV,VEC>& a)
{typename TV::SCALAR t=tan(a.x),s=1+t*t;return Make_Diff(t,s*a.dx);}

template<class TV,class VEC> inline AUTO_DIFF_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC())))>
abs(const AUTO_DIFF_EXT<TV,VEC>& a)
{
    AUTO_DIFF_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC())))> r;
    if(a.x>=0) r.Fill_From(a);
    else r.Fill_From(-a);
    return r;
}

template<class TV,class VEC,class VEC1> inline
decltype(Make_Diff(SC,SC*DX+SC*DX1))
hypot(const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC1>& b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return Make_Diff(c,d*a.dx+e*b.dx);
}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
hypot(const AUTO_DIFF_EXT<TV,VEC>& a,typename TV::SCALAR b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return Make_Diff(c,d*a.dx);
}

template<class TV,class VEC> inline decltype(hypot(AUTO_DIFF_EXT<TV,VEC>(),SC))
hypot(typename TV::SCALAR b,const AUTO_DIFF_EXT<TV,VEC>& a)
{return hypot(a,b);}

template<class TV,class VEC,class VEC1,class VEC2> inline
decltype(Make_Diff(SC,SC*DX+SC*DX1+SC*DX1))
hypot(const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC1>& b,const AUTO_DIFF_EXT<TV,VEC2>& c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(bb*c.dx-cc*b.dx);
    auto ca=Outer_Product(cc*a.dx-aa*c.dx);
    return Make_Diff(s,aa*a.dx+bb*b.dx+cc*c.dx);
}

template<class TV,class VEC,class VEC1> inline
decltype(Make_Diff(SC,SC*DX+SC*DX1))
hypot(const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC1>& b,typename TV::SCALAR c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(b.dx);
    auto ca=Outer_Product(a.dx);
    return Make_Diff(s,aa*a.dx+bb*b.dx);
}

template<class TV,class VEC,class VEC1> inline decltype(hypot(AUTO_DIFF_EXT<TV,VEC>(),AUTO_DIFF_EXT<TV,VEC1>(),SC))
hypot(const AUTO_DIFF_EXT<TV,VEC>& a,typename TV::SCALAR c,const AUTO_DIFF_EXT<TV,VEC1>& b)
{return hypot(a,b,c);}

template<class TV,class VEC,class VEC1> inline decltype(hypot(AUTO_DIFF_EXT<TV,VEC>(),AUTO_DIFF_EXT<TV,VEC1>(),SC))
hypot(typename TV::SCALAR c,const AUTO_DIFF_EXT<TV,VEC>& a,const AUTO_DIFF_EXT<TV,VEC1>& b)
{return hypot(a,b,c);}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
hypot(const AUTO_DIFF_EXT<TV,VEC>& a,typename TV::SCALAR b,typename TV::SCALAR c)
{
    typename TV::SCALAR t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return Make_Diff(s,d*a.dx);
}

template<class TV,class VEC> inline decltype(hypot(AUTO_DIFF_EXT<TV,VEC>(),SC,SC))
hypot(typename TV::SCALAR b,const AUTO_DIFF_EXT<TV,VEC>& a,typename TV::SCALAR c)
{return hypot(a,b,c);}

template<class TV,class VEC> inline decltype(hypot(AUTO_DIFF_EXT<TV,VEC>(),SC,SC))
hypot(typename TV::SCALAR b,typename TV::SCALAR c,const AUTO_DIFF_EXT<TV,VEC>& a)
{return hypot(a,b,c);}

template<class TV,class VEC,class VEC1> inline decltype(Make_Diff(SC,SC*DX-SC*DX1))
atan2(const AUTO_DIFF_EXT<TV,VEC>& y,const AUTO_DIFF_EXT<TV,VEC1>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    auto f=d*y.dx-e*x.dx;
    return Make_Diff(::std::atan2(y.x,x.x),f);
}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
atan2(const AUTO_DIFF_EXT<TV,VEC>& y,typename TV::SCALAR x)
{
    typename TV::SCALAR c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    auto f=d*y.dx;
    return Make_Diff(::std::atan2(y.x,x),f);
}

template<class TV,class VEC> inline decltype(Make_Diff(SC,SC*DX))
atan2(typename TV::SCALAR y,const AUTO_DIFF_EXT<TV,VEC>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    auto f=-e*x.dx;
    return Make_Diff(::std::atan2(y,x.x),f);
}
template<class TV,class VEC> struct AUTO_DIFF_EXT_VEC;

template<class TV,class VEC> AUTO_DIFF_EXT_VEC<TV,VEC>
Make_Diff_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx);

template<class TV,class VEC> struct AUTO_DIFF_EXT_VEC;

template<class T,class TV,class VEC1,class VEC2>
decltype(Make_Diff_Vec(TV(),MATRIX<T,3>()*VDX2-MATRIX<T,3>()*VDX1))
Cross_Helper(const AUTO_DIFF_EXT_VEC<VECTOR<T,3>,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC2>& b);

template<class TV,class VEC1,class VEC2>
typename DISABLE_IF<TV::m==3,const AUTO_DIFF_EXT_VEC<TV,VEC1> >::TYPE
Cross_Helper(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC2>& b);

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
    decltype(Make_Diff_Vec(TV(),dx+VDX1)) operator+(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const
    {return Make_Diff_Vec(x+a.x,dx+a.dx);}

    template<class VEC1>
    decltype(Make_Diff_Vec(TV(),dx-VDX1)) operator-(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const
    {return Make_Diff_Vec(x-a.x,dx-a.dx);}

    template<class VEC1>
    decltype(Make_Diff_Vec(x*SC,SC*dx+Outer_Product(x,DX1)))
    operator*(const AUTO_DIFF_EXT<TV,VEC1>& a) const
    {return Make_Diff_Vec(x*a.x,a.x*dx+Outer_Product(x,a.dx));}

    template<class VEC1>
    decltype(Make_Diff_Vec(x/SC,dx/SC+Outer_Product(x,-DX1/(SC*SC))))
    operator/(const AUTO_DIFF_EXT<TV,VEC1>& a) const
    {auto p=-a.dx/(a.x*a.x);return Make_Diff_Vec(x/a.x,dx/a.x+Outer_Product(x,p));}

    AUTO_DIFF_EXT_VEC operator+(TV a) const
    {return Make_Diff_Vec(x+a,dx);}

    AUTO_DIFF_EXT_VEC operator-(TV a) const
    {return Make_Diff_Vec(x-a,dx);}

    decltype(Make_Diff_Vec(x*SC,dx*SC)) operator*(T a) const
    {return Make_Diff_Vec(x*a,dx*a);}

    decltype(Make_Diff_Vec(x/SC,dx/SC)) operator/(T a) const
    {return Make_Diff_Vec(x/a,dx/a);}

    decltype(Make_Diff(x.Dot(TV()),dx.Transpose_Times(TV()))) Dot(TV v) const
    {return Make_Diff(x.Dot(v),dx.Transpose_Times(v),Contract_0);}

    template<class VEC1>
    decltype(Make_Diff(x.Dot(TV()),dx.Transpose_Times(TV())+VDX1.Transpose_Times(x)))
    Dot(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const
    {return Make_Diff(x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x));}

    decltype(Make_Diff_Vec(TV(),-MATRIX<T,TV::m>()*dx))
    Cross(const TV& a) const
    {
        MATRIX<T,TV::m> cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(a);
        return Make_Diff_Vec(x.Cross(a),-cp_a*dx,Contract_0);
    }

    template<class VEC1>
    auto Cross(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a) const -> decltype(Cross_Helper(*this,AUTO_DIFF_EXT_VEC<TV,VEC1>()))
    {return Cross_Helper(*this,a);}

    decltype(Make_Diff(x.Magnitude_Squared(),dx.Transpose_Times(x)*2)) Magnitude_Squared() const
    {return Make_Diff(x.Magnitude_Squared(),dx.Transpose_Times(x)*2);}

    decltype(Make_Diff(SC,dx.Transpose_Times(x)/SC)) Magnitude() const
    {T s=x.Magnitude();auto t=dx.Transpose_Times(x)/s;return Make_Diff(s,t);}
};


template<class T,class TV,class VEC1,class VEC2>
decltype(Make_Diff_Vec(TV(),MATRIX<T,3>()*VDX2-MATRIX<T,3>()*VDX1))
Cross_Helper(const AUTO_DIFF_EXT_VEC<VECTOR<T,3>,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC2>& b)
{
    MATRIX<T,TV::m> cp_t=MATRIX<T,TV::m>::Cross_Product_Matrix(a.x),cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(b.x);
    return Make_Diff_Vec(a.x.Cross(b.x),cp_t*b.dx-cp_a*a.dx);
}

template<class TV,class VEC1,class VEC2>
typename DISABLE_IF<TV::m==3,const AUTO_DIFF_EXT_VEC<TV,VEC1> >::TYPE
Cross_Helper(const AUTO_DIFF_EXT_VEC<TV,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC2>& b)
{
    PHYSBAM_FATAL_ERROR("Cross product not defined except in 3D.");
    return a;
}

template<class TV,class VEC> inline AUTO_DIFF_EXT_VEC<TV,VEC>
Make_Diff_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx)
{return AUTO_DIFF_EXT_VEC<TV,VEC>(x,dx);}

template<int n,int i,class TV> inline
AUTO_DIFF_EXT_VEC<TV,typename ONE_NONZERO_VECTOR_MAT<TV,n,i>::TYPE> Diff_From_Var(TV v)
{return AUTO_DIFF_EXT_VEC<TV,typename ONE_NONZERO_VECTOR_MAT<TV,n,i>::TYPE>(v);}

template<class TV,class VEC> inline decltype(Make_Diff_Vec(TV(),Outer_Product(TV(),DX)))
operator*(const AUTO_DIFF_EXT<TV,VEC>& a,TV v)
{return Make_Diff_Vec(v*a.x,Outer_Product(v,a.dx),Tensor_Product_0);}

template<class TV,class VEC> inline decltype(AUTO_DIFF_EXT<TV,VEC>()*TV())
operator*(TV v,const AUTO_DIFF_EXT<TV,VEC>& a)
{return a*v;}

template<class TV,class VEC,class VEC1> inline
decltype(AUTO_DIFF_EXT_VEC<TV,VEC>()*AUTO_DIFF_EXT<TV,VEC1>())
operator*(const AUTO_DIFF_EXT<TV,VEC1>& a,const AUTO_DIFF_EXT_VEC<TV,VEC>& v)
{return v*a;}

template<class TV,class VEC>
decltype(AUTO_DIFF_EXT_VEC<TV,VEC>()*SC)
operator*(typename TV::SCALAR a,const AUTO_DIFF_EXT_VEC<TV,VEC>& v)
{return v*a;}

template<class TV,class VEC> AUTO_DIFF_EXT_VEC<TV,VEC>
operator+(TV a,const AUTO_DIFF_EXT_VEC<TV,VEC>& u)
{return Make_Diff_Vec(u.x+a,u.dx);}

template<class TV,class VEC> decltype(Make_Diff_Vec(TV(),-VDX))
operator-(TV a,const AUTO_DIFF_EXT_VEC<TV,VEC>& u)
{return Make_Diff_Vec(a-u.x,-u.dx);}

template<class TV,class VEC1>
decltype(TV()*(SC/AUTO_DIFF_EXT<TV,VEC1>()))
operator/(TV u,const AUTO_DIFF_EXT<TV,VEC1>& a)
{return u*((typename TV::SCALAR)1/a);}

template<class TV,class VEC>
decltype(Make_Diff_Vec(MATRIX<typename TV::SCALAR,TV::m>()*TV(),MATRIX<typename TV::SCALAR,TV::m>()*VDX))
operator*(const MATRIX<typename TV::SCALAR,TV::m>& m,const AUTO_DIFF_EXT_VEC<TV,VEC>& v)
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
