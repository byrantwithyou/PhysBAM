//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_HESS_EXT
//##################################################################### 
#ifndef __AUTO_HESS_EXT__
#define __AUTO_HESS_EXT__

#include <Tools/Auto_Diff/HESSIAN.h>
#include <Tools/Auto_Diff/HESSIAN_VEC.h>
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

#define SC typename TV::SCALAR()
#define DX GRADIENT<TV,VEC>()
#define DDX HESSIAN<TV,MAT>()
#define DX1 GRADIENT<TV,VEC1>()
#define DDX1 HESSIAN<TV,MAT1>()
#define DX2 GRADIENT<TV,VEC2>()
#define DDX2 HESSIAN<TV,MAT2>()
#define VDX GRADIENT_VEC<TV,VEC>()
#define VDDX HESSIAN_VEC<TV,MAT>()
#define VDX1 GRADIENT_VEC<TV,VEC1>()
#define VDDX1 HESSIAN_VEC<TV,MAT1>()
#define VDX2 GRADIENT_VEC<TV,VEC2>()
#define VDDX2 HESSIAN_VEC<TV,MAT2>()
template<class TV,class VEC,class MAT> struct AUTO_HESS_EXT;

template<class TV,class VEC,class MAT> AUTO_HESS_EXT<TV,VEC,MAT>
Make_Hess(typename TV::SCALAR x,const GRADIENT<TV,VEC>& dx,const HESSIAN<TV,MAT>& ddx);

template<class TV,class VEC,class MAT>
struct AUTO_HESS_EXT
{
    typedef typename TV::SCALAR T;

    T x;
    GRADIENT<TV,VEC> dx;
    HESSIAN<TV,MAT> ddx;

    AUTO_HESS_EXT(T z=T()): x(z)
    {}

    AUTO_HESS_EXT(T z,const GRADIENT<TV,VEC>& dz,const HESSIAN<TV,MAT>& ddz)
        :x(z),dx(dz),ddx(ddz)
    {}

    decltype(Make_Hess(-x,-dx,-ddx)) operator-() const
    {return Make_Hess(-x,-dx,-ddx);}

    AUTO_HESS_EXT operator+() const
    {return *this;}

    template<class VEC1,class MAT1> decltype(Make_Hess(x+SC,dx+DX1,ddx+DDX1))
    operator+(const AUTO_HESS_EXT<TV,VEC1,MAT1>& a) const
    {return Make_Hess(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class VEC1,class MAT1> decltype(Make_Hess(x-SC,dx-DX1,ddx-DDX1))
    operator-(const AUTO_HESS_EXT<TV,VEC1,MAT1>& a) const
    {return Make_Hess(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class VEC1,class MAT1> decltype(Make_Hess(x*SC,SC*dx+x*DX1,SC*ddx+x*DDX1+Symmetric_Outer_Product(dx,DX1)))
    operator*(const AUTO_HESS_EXT<TV,VEC1,MAT1>& a) const
    {return Make_Hess(x*a.x,a.x*dx+x*a.dx,a.x*ddx+x*a.ddx+Symmetric_Outer_Product(dx,a.dx));}

    template<class VEC1,class MAT1> decltype(Make_Hess(x/SC,(dx-(x/SC)*DX1)/SC,ddx/SC-x/SC*DDX1/SC-Symmetric_Outer_Product((dx-(x/SC)*DX1),DX1/(SC*SC))))
    operator/(const AUTO_HESS_EXT<TV,VEC1,MAT1>& a) const
    {T z=x/a.x,ax2=sqr(a.x);auto q=dx-z*a.dx;return Make_Hess(z,q/a.x,ddx/a.x-z*a.ddx/a.x-Symmetric_Outer_Product(q,a.dx/ax2));}

    AUTO_HESS_EXT operator+(T a) const
    {return Make_Hess(x+a,dx,ddx);}

    AUTO_HESS_EXT operator-(T a) const
    {return Make_Hess(x-a,dx,ddx);}

    decltype(Make_Hess(x*SC,dx*SC,ddx*SC)) operator*(T a) const
    {return Make_Hess(x*a,a*dx,a*ddx);}

    decltype(Make_Hess(x/SC,dx/SC,ddx/SC)) operator/(T a) const
    {return Make_Hess(x/a,dx/a,ddx/a);}

    template<class VEC1,class MAT1>
    void Fill_From(const AUTO_HESS_EXT<TV,VEC1,MAT1>& z)
    {x=z.x;::PhysBAM::HETERO_DIFF::Fill_From(dx,z.dx);::PhysBAM::HETERO_DIFF::Fill_From(ddx,z.ddx);}
};

template<class TV,class VEC,class MAT> AUTO_HESS_EXT<TV,VEC,MAT>
Make_Hess(typename TV::SCALAR x,const GRADIENT<TV,VEC>& dx,const HESSIAN<TV,MAT>& ddx)
{return AUTO_HESS_EXT<TV,VEC,MAT>(x,dx,ddx);}

template<int n,int i,class TV>
AUTO_HESS_EXT<TV,typename ONE_NONZERO_VECTOR<TV,n,i>::TYPE,typename EMPTY_MAT<TV,n>::TYPE>
From_Var(TV v,int j)
{
    AUTO_HESS_EXT<TV,typename ONE_NONZERO_VECTOR<TV,n,i>::TYPE,typename EMPTY_MAT<TV,n>::TYPE> r(v(j));
    Set<i>(r.dx.x,TV::Axis_Vector(j));
    return r;
}

template<class TV,int n>
AUTO_HESS_EXT<TV,typename EMPTY_VEC<TV,n>::TYPE,typename EMPTY_MAT<TV,n>::TYPE>
From_Const(typename TV::SCALAR a)
{return AUTO_HESS_EXT<TV,typename EMPTY_VEC<TV,n>::TYPE,typename EMPTY_MAT<TV,n>::TYPE>(a);}

template<class TV,class VEC,class MAT>
inline AUTO_HESS_EXT<TV,VEC,MAT> operator+(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{return Make_Hess(c+a.x,a.dx,a.ddx);}

template<class TV,class VEC,class MAT>
inline AUTO_HESS_EXT<TV,VEC,MAT> operator-(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{return Make_Hess(c-a.x,-a.dx,-a.ddx);}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*DDX))
operator*(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{return Make_Hess(c*a.x,c*a.dx,c*a.ddx);}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*Outer_Product(DX)-SC*DDX))
operator/(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{typename TV::SCALAR z=c/a.x,w=z/a.x;return Make_Hess(z,-w*a.dx,2*w/a.x*Outer_Product(a.dx)-w*a.ddx);}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,DX/SC,DDX/SC-Outer_Product(DX/SC)/SC))
sqrt(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{typename TV::SCALAR s=sqrt(a.x);auto t=a.dx/(2*s);return Make_Hess(s,t,a.ddx/(2*s)-Outer_Product(t)/s);}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(sqr(SC),2*SC*DX,2*SC*DDX+Outer_Product(DX)*2))
sqr(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{return Make_Hess(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+Outer_Product(a.dx)*2);}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(cube(SC),SC*DX,SC*DDX+Outer_Product(DX)*SC))
cube(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{typename TV::SCALAR sq=sqr(a.x);return Make_Hess(a.x*sq,3*sq*a.dx,3*sq*a.ddx+Outer_Product(a.dx)*(6*a.x));}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline AUTO_HESS_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))>
max(const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC1,MAT1>& b)
{
    AUTO_HESS_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))> r;
    if(a.x>b.x) r.Fill_From(a);
    else r.Fill_From(b);
    return r;
}

template<class TV,class VEC,class MAT> inline decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC)
max(const AUTO_HESS_EXT<TV,VEC,MAT>& a,typename TV::SCALAR b)
{
    decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class TV,class VEC,class MAT> inline decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC)
max(typename TV::SCALAR b,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{
    decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class TV,class VEC,class MAT>
inline AUTO_HESS_EXT<TV,VEC,MAT>
max(const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC,MAT>& b)
{return a.x>b.x?a:b;}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline AUTO_HESS_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))>
min(const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC1,MAT1>& b)
{
    AUTO_HESS_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))> r;
    if(a.x>b.x) r.Fill_From(b);
    else r.Fill_From(a);
    return r;
}

template<class TV,class VEC,class MAT> inline decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC)
min(const AUTO_HESS_EXT<TV,VEC,MAT>& a,typename TV::SCALAR b)
{
    decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class TV,class VEC,class MAT> inline decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC)
min(typename TV::SCALAR b,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{
    decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*SC) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class TV,class VEC,class MAT>
inline AUTO_HESS_EXT<TV,VEC,MAT>
min(const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC,MAT>& b)
{return a.x>b.x?b:a;}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,DX/SC,DDX/SC-Outer_Product(DX/SC)))
log(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{auto z=a.dx/a.x;return Make_Hess(::std::log(a.x),z,a.ddx/a.x-Outer_Product(z));}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*(DDX+Outer_Product(DX))))
exp(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{typename TV::SCALAR s=exp(a.x);return Make_Hess(s,s*a.dx,s*(a.ddx+Outer_Product(a.dx)));}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*DDX-SC*Outer_Product(DX)))
sin(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return Make_Hess(s,c*a.dx,c*a.ddx-s*Outer_Product(a.dx));}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,-SC*DX,-SC*DDX-SC*Outer_Product(DX)))
cos(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return Make_Hess(c,-s*a.dx,-s*a.ddx-c*Outer_Product(a.dx));}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*DDX+SC*Outer_Product(DX)))
tan(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{typename TV::SCALAR t=tan(a.x),s=1+t*t;return Make_Hess(t,s*a.dx,s*a.ddx+2*s*t*Outer_Product(a.dx));}

template<class TV,class VEC,class MAT> inline AUTO_HESS_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC()))),decltype(MAT_CHOOSE::Type(MAT(),MAT_NEG::Type(MAT())))>
abs(const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{
    AUTO_HESS_EXT<TV,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC()))),decltype(MAT_CHOOSE::Type(MAT(),MAT_NEG::Type(MAT())))> r;
    if(a.x>=0) r.Fill_From(a);
    else r.Fill_From(-a);
    return r;
}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline
decltype(Make_Hess(SC,SC*DX+SC*DX1,SC*DDX+SC*DDX1+Outer_Product(SC*DX1-SC*DX)/SC))
hypot(const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC1,MAT1>& b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return Make_Hess(c,d*a.dx+e*b.dx,d*a.ddx+e*b.ddx+Outer_Product(d*b.dx-e*a.dx)/c);
}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*DDX+Outer_Product(SC*DX)/SC))
hypot(const AUTO_HESS_EXT<TV,VEC,MAT>& a,typename TV::SCALAR b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return Make_Hess(c,d*a.dx,d*a.ddx+Outer_Product(e*a.dx)/c);
}

template<class TV,class VEC,class MAT> inline decltype(hypot(AUTO_HESS_EXT<TV,VEC,MAT>(),SC))
hypot(typename TV::SCALAR b,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{return hypot(a,b);}

template<class TV,class VEC,class MAT,class VEC1,class MAT1,class VEC2,class MAT2> inline
decltype(Make_Hess(SC,SC*DX+SC*DX1+SC*DX1,SC*DDX+SC*DDX1+SC*DDX1+(Outer_Product(SC*DX1-SC*DX)+Outer_Product(SC*DX1-SC*DX1)+Outer_Product(SC*DX-SC*DX1))/SC))
hypot(const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC1,MAT1>& b,const AUTO_HESS_EXT<TV,VEC2,MAT2>& c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(bb*c.dx-cc*b.dx);
    auto ca=Outer_Product(cc*a.dx-aa*c.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx+cc*c.dx,aa*a.ddx+bb*b.ddx+cc*c.ddx+(ab+bc+ca)/s);
}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline
decltype(Make_Hess(SC,SC*DX+SC*DX1,SC*DDX+SC*DDX1+Outer_Product(SC*DX1-SC*DX)/SC+(Outer_Product(DX1)+Outer_Product(DX))*SC))
hypot(const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC1,MAT1>& b,typename TV::SCALAR c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(b.dx);
    auto ca=Outer_Product(a.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx,aa*a.ddx+bb*b.ddx+ab/s+(bc+ca)*(sqr(cc)/s));
}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline decltype(hypot(AUTO_HESS_EXT<TV,VEC,MAT>(),AUTO_HESS_EXT<TV,VEC1,MAT1>(),SC))
hypot(const AUTO_HESS_EXT<TV,VEC,MAT>& a,typename TV::SCALAR c,const AUTO_HESS_EXT<TV,VEC1,MAT1>& b)
{return hypot(a,b,c);}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline decltype(hypot(AUTO_HESS_EXT<TV,VEC,MAT>(),AUTO_HESS_EXT<TV,VEC1,MAT1>(),SC))
hypot(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,VEC,MAT>& a,const AUTO_HESS_EXT<TV,VEC1,MAT1>& b)
{return hypot(a,b,c);}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*DDX+Outer_Product(DX)*SC))
hypot(const AUTO_HESS_EXT<TV,VEC,MAT>& a,typename TV::SCALAR b,typename TV::SCALAR c)
{
    typename TV::SCALAR t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return Make_Hess(s,d*a.dx,d*a.ddx+Outer_Product(a.dx)*e);
}

template<class TV,class VEC,class MAT> inline decltype(hypot(AUTO_HESS_EXT<TV,VEC,MAT>(),SC,SC))
hypot(typename TV::SCALAR b,const AUTO_HESS_EXT<TV,VEC,MAT>& a,typename TV::SCALAR c)
{return hypot(a,b,c);}

template<class TV,class VEC,class MAT> inline decltype(hypot(AUTO_HESS_EXT<TV,VEC,MAT>(),SC,SC))
hypot(typename TV::SCALAR b,typename TV::SCALAR c,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{return hypot(a,b,c);}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline decltype(Make_Hess(SC,SC*DX-SC*DX1,SC*DDX-SC*DDX1-Symmetric_Outer_Product(SC*DX-SC*DX1,SC*DX1+SC*DX)))
atan2(const AUTO_HESS_EXT<TV,VEC,MAT>& y,const AUTO_HESS_EXT<TV,VEC1,MAT1>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    auto f=d*y.dx-e*x.dx;
    return Make_Hess(::std::atan2(y.x,x.x),f,d*y.ddx-e*x.ddx-Symmetric_Outer_Product(f,d*x.dx+e*y.dx));
}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*DDX-Symmetric_Outer_Product(SC*DX,SC*DX)))
atan2(const AUTO_HESS_EXT<TV,VEC,MAT>& y,typename TV::SCALAR x)
{
    typename TV::SCALAR c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    auto f=d*y.dx;
    return Make_Hess(::std::atan2(y.x,x),f,d*y.ddx-Symmetric_Outer_Product(f,e*y.dx));
}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess(SC,SC*DX,SC*DDX-Symmetric_Outer_Product(SC*DX,SC*DX)))
atan2(typename TV::SCALAR y,const AUTO_HESS_EXT<TV,VEC,MAT>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    auto f=-e*x.dx;
    return Make_Hess(::std::atan2(y,x.x),f,-e*x.ddx-Symmetric_Outer_Product(f,d*x.dx));
}
template<class TV,class VEC,class MAT> struct AUTO_HESS_EXT_VEC;

template<class TV,class VEC,class MAT> AUTO_HESS_EXT_VEC<TV,VEC,MAT>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const HESSIAN_VEC<TV,MAT>& ddx);

template<class TV,class VEC,class MAT> struct AUTO_HESS_EXT_VEC;

template<class T,class TV,class VEC1,class MAT1,class VEC2,class MAT2>
decltype(Make_Hess_Vec(TV(),MATRIX<T,3>()*VDX2-MATRIX<T,3>()*VDX1,Contract_0(VDDX1,MATRIX<T,3>())-Contract_0(VDDX2,MATRIX<T,3>())+Symmetric_Double_Contract_12_With_Tensor(PERMUTATION_TENSOR<T>(1),VDX1,VDX2)))
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2>& b);

template<class TV,class VEC1,class MAT1,class VEC2,class MAT2>
typename DISABLE_IF<TV::m==3,const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1> >::TYPE
Cross_Helper(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2>& b);

template<class TV,class VEC,class MAT>
struct AUTO_HESS_EXT_VEC
{
    typedef typename TV::SCALAR T;

    TV x;
    GRADIENT_VEC<TV,VEC> dx;
    HESSIAN_VEC<TV,MAT> ddx;

    AUTO_HESS_EXT_VEC(TV x=TV()):
        x(x)
    {}

    AUTO_HESS_EXT_VEC(TV x,const GRADIENT_VEC<TV,VEC>& dx,const HESSIAN_VEC<TV,MAT>& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    decltype(Make_Hess_Vec(-TV(),-GRADIENT_VEC<TV,VEC>(),-HESSIAN_VEC<TV,MAT>())) operator-() const
    {return Make_Hess_Vec(-x,-dx,-ddx);}

    AUTO_HESS_EXT_VEC operator+() const
    {return *this;}

    template<class VEC1,class MAT1>
    decltype(Make_Hess_Vec(TV(),dx+VDX1,ddx+VDDX1)) operator+(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const
    {return Make_Hess_Vec(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class VEC1,class MAT1>
    decltype(Make_Hess_Vec(TV(),dx-VDX1,ddx-VDDX1)) operator-(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const
    {return Make_Hess_Vec(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class VEC1,class MAT1>
    decltype(Make_Hess_Vec(x*SC,SC*dx+Outer_Product(x,DX1),ddx*SC+Tensor_Product_0(DDX1,x)+Symmetric_Tensor_Product_12(dx,DX1)))
    operator*(const AUTO_HESS_EXT<TV,VEC1,MAT1>& a) const
    {return Make_Hess_Vec(x*a.x,a.x*dx+Outer_Product(x,a.dx),ddx*a.x+Tensor_Product_0(a.ddx,x)+Symmetric_Tensor_Product_12(dx,a.dx));}

    template<class VEC1,class MAT1>
    decltype(Make_Hess_Vec(x/SC,dx/SC+Outer_Product(x,-DX1/(SC*SC)),ddx/SC+Tensor_Product_0(2/SC*Outer_Product(DX1)-DDX1,x/(SC*SC))+Symmetric_Tensor_Product_12(dx,-DX1/(SC*SC))))
    operator/(const AUTO_HESS_EXT<TV,VEC1,MAT1>& a) const
    {auto p=-a.dx/(a.x*a.x);return Make_Hess_Vec(x/a.x,dx/a.x+Outer_Product(x,p),ddx/a.x+Tensor_Product_0(2/a.x*Outer_Product(a.dx)-a.ddx,x/(a.x*a.x))+Symmetric_Tensor_Product_12(dx,p));}

    AUTO_HESS_EXT_VEC operator+(TV a) const
    {return Make_Hess_Vec(x+a,dx,ddx);}

    AUTO_HESS_EXT_VEC operator-(TV a) const
    {return Make_Hess_Vec(x-a,dx,ddx);}

    decltype(Make_Hess_Vec(x*SC,dx*SC,ddx*SC)) operator*(T a) const
    {return Make_Hess_Vec(x*a,dx*a,ddx*a);}

    decltype(Make_Hess_Vec(x/SC,dx/SC,ddx/SC)) operator/(T a) const
    {return Make_Hess_Vec(x/a,dx/a,ddx/a);}

    decltype(Make_Hess(x.Dot(TV()),dx.Transpose_Times(TV()),Contract_0(ddx,TV()))) Dot(TV v) const
    {return Make_Hess(x.Dot(v),dx.Transpose_Times(v),Contract_0(ddx,v));}

    template<class VEC1,class MAT1>
    decltype(Make_Hess(x.Dot(TV()),dx.Transpose_Times(TV())+VDX1.Transpose_Times(x),Contract_0(ddx,TV())+Contract_0(VDDX1,x)+Symmetric_Transpose_Times(dx,VDX1)))
    Dot(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const
    {return Make_Hess(x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x),Contract_0(ddx,a.x)+Contract_0(a.ddx,x)+Symmetric_Transpose_Times(dx,a.dx));}

    decltype(Make_Hess_Vec(TV(),-MATRIX<T,TV::m>()*dx,Contract_0(ddx,MATRIX<T,TV::m>())))
    Cross(const TV& a) const
    {
        MATRIX<T,TV::m> cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(a);
        return Make_Hess_Vec(x.Cross(a),-cp_a*dx,Contract_0(ddx,cp_a));
    }

    template<class VEC1,class MAT1>
    auto Cross(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const -> decltype(Cross_Helper(*this,AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>()))
    {return Cross_Helper(*this,a);}

    decltype(Make_Hess(x.Magnitude_Squared(),dx.Transpose_Times(x)*2,(Contract_0(ddx,x)+Transpose_Times_Self(dx))*2)) Magnitude_Squared() const
    {return Make_Hess(x.Magnitude_Squared(),dx.Transpose_Times(x)*2,(Contract_0(ddx,x)+Transpose_Times_Self(dx))*2);}

    decltype(Make_Hess(SC,dx.Transpose_Times(x)/SC,(Contract_0(ddx,x)+Transpose_Times_Self(dx))/SC-Outer_Product(dx.Transpose_Times(x)/SC)/SC)) Magnitude() const
    {T s=x.Magnitude();auto t=dx.Transpose_Times(x)/s;return Make_Hess(s,t,(Contract_0(ddx,x)+Transpose_Times_Self(dx))/s-Outer_Product(t)/s);}
};


template<class T,class TV,class VEC1,class MAT1,class VEC2,class MAT2>
decltype(Make_Hess_Vec(TV(),MATRIX<T,3>()*VDX2-MATRIX<T,3>()*VDX1,Contract_0(VDDX1,MATRIX<T,3>())-Contract_0(VDDX2,MATRIX<T,3>())+Symmetric_Double_Contract_12_With_Tensor(PERMUTATION_TENSOR<T>(1),VDX1,VDX2)))
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2>& b)
{
    MATRIX<T,TV::m> cp_t=MATRIX<T,TV::m>::Cross_Product_Matrix(a.x),cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(b.x);
    return Make_Hess_Vec(a.x.Cross(b.x),cp_t*b.dx-cp_a*a.dx,Contract_0(a.ddx,cp_a)-Contract_0(b.ddx,cp_t)+Symmetric_Double_Contract_12_With_Tensor(PERMUTATION_TENSOR<T>(1),a.dx,b.dx));
}

template<class TV,class VEC1,class MAT1,class VEC2,class MAT2>
typename DISABLE_IF<TV::m==3,const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1> >::TYPE
Cross_Helper(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2>& b)
{
    PHYSBAM_FATAL_ERROR("Cross product not defined except in 3D.");
    return a;
}

template<class TV,class VEC,class MAT> inline AUTO_HESS_EXT_VEC<TV,VEC,MAT>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const HESSIAN_VEC<TV,MAT>& ddx)
{return AUTO_HESS_EXT_VEC<TV,VEC,MAT>(x,dx,ddx);}

template<int n,int i,class TV> inline
AUTO_HESS_EXT_VEC<TV,typename ONE_NONZERO_VECTOR_MAT<TV,n,i>::TYPE,typename EMPTY_MAT_TEN<TV,n>::TYPE> From_Var(TV v)
{return AUTO_HESS_EXT_VEC<TV,typename ONE_NONZERO_VECTOR_MAT<TV,n,i>::TYPE,typename EMPTY_MAT_TEN<TV,n>::TYPE>(v);}

template<class TV,class VEC,class MAT> inline decltype(Make_Hess_Vec(TV(),Outer_Product(TV(),DX),Tensor_Product_0(DDX,TV())))
operator*(const AUTO_HESS_EXT<TV,VEC,MAT>& a,TV v)
{return Make_Hess_Vec(v*a.x,Outer_Product(v,a.dx),Tensor_Product_0(a.ddx,v));}

template<class TV,class VEC,class MAT> inline decltype(AUTO_HESS_EXT<TV,VEC,MAT>()*TV())
operator*(TV v,const AUTO_HESS_EXT<TV,VEC,MAT>& a)
{return a*v;}

template<class TV,class VEC,class MAT,class VEC1,class MAT1> inline
decltype(AUTO_HESS_EXT_VEC<TV,VEC,MAT>()*AUTO_HESS_EXT<TV,VEC1,MAT1>())
operator*(const AUTO_HESS_EXT<TV,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& v)
{return v*a;}

template<class TV,class VEC,class MAT>
decltype(AUTO_HESS_EXT_VEC<TV,VEC,MAT>()*SC)
operator*(typename TV::SCALAR a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& v)
{return v*a;}

template<class TV,class VEC,class MAT> AUTO_HESS_EXT_VEC<TV,VEC,MAT>
operator+(TV a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& u)
{return Make_Hess_Vec(u.x+a,u.dx,u.ddx);}

template<class TV,class VEC,class MAT> decltype(Make_Hess_Vec(TV(),-VDX,-VDDX))
operator-(TV a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& u)
{return Make_Hess_Vec(a-u.x,-u.dx,-u.ddx);}

template<class TV,class VEC1,class MAT1>
decltype(TV()*(SC/AUTO_HESS_EXT<TV,VEC1,MAT1>()))
operator/(TV u,const AUTO_HESS_EXT<TV,VEC1,MAT1>& a)
{return u*((typename TV::SCALAR)1/a);}

template<class TV,class VEC,class MAT>
decltype(Make_Hess_Vec(MATRIX<typename TV::SCALAR,TV::m>()*TV(),MATRIX<typename TV::SCALAR,TV::m>()*VDX,Contract_0(VDDX,MATRIX<typename TV::SCALAR,TV::m>())))
operator*(const MATRIX<typename TV::SCALAR,TV::m>& m,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& v)
{return Make_Hess_Vec(m*v.x,m*v.dx,Contract_0(v.ddx,m.Transposed()));}
}
using HETERO_DIFF::AUTO_HESS_EXT;
using HETERO_DIFF::AUTO_HESS_EXT_VEC;
using HETERO_DIFF::From_Const;
using HETERO_DIFF::From_Var;
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
