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
using ::PhysBAM::cube;

template<class T,class VEC,class MAT> struct AUTO_HESS_EXT;

template<class T,class VEC,class MAT> AUTO_HESS_EXT<T,VEC,MAT>
Make_Hess(T x,const GRADIENT<T,VEC>& dx,const HESSIAN<T,MAT>& ddx);

template<class T,class VEC,class MAT>
struct AUTO_HESS_EXT
{
    static_assert(is_scalar<T>::value,"This definition of AUTO_HESS_EXT is only for scalars");

    T x;
    GRADIENT<T,VEC> dx;
    HESSIAN<T,MAT> ddx;

    AUTO_HESS_EXT(T z=T()): x(z)
    {}

    AUTO_HESS_EXT(T z,const GRADIENT<T,VEC>& dz,const HESSIAN<T,MAT>& ddz)
        :x(z),dx(dz),ddx(ddz)
    {}

    decltype(Make_Hess(-x,-dx,-ddx)) operator-() const
    {return Make_Hess(-x,-dx,-ddx);}

    AUTO_HESS_EXT operator+() const
    {return *this;}

    template<class VEC1,class MAT1> auto
    operator+(const AUTO_HESS_EXT<T,VEC1,MAT1>& a) const
        -> decltype(Make_Hess(this->x+a.x,this->dx+a.dx,this->ddx+a.ddx))
    {return Make_Hess(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class VEC1,class MAT1> auto
    operator-(const AUTO_HESS_EXT<T,VEC1,MAT1>& a) const
        -> decltype(Make_Hess(x-a.x,dx-a.dx,ddx-a.ddx))
    {return Make_Hess(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class VEC1,class MAT1> auto
    operator*(const AUTO_HESS_EXT<T,VEC1,MAT1>& a) const
        -> decltype(Make_Hess(this->x*a.x,a.x*this->dx+this->x*a.dx,a.x*this->ddx+this->x*a.ddx+Symmetric_Outer_Product(this->dx,a.dx)))
    {return Make_Hess(x*a.x,a.x*dx+x*a.dx,a.x*ddx+x*a.ddx+Symmetric_Outer_Product(dx,a.dx));}

    template<class VEC1,class MAT1> auto
    operator/(const AUTO_HESS_EXT<T,VEC1,MAT1>& a) const
        -> decltype(Make_Hess(T(),(this->dx-T()*a.dx)/a.x,this->ddx/a.x-T()*a.ddx/a.x-Symmetric_Outer_Product((this->dx-T()*a.dx),a.dx/T())))
    {
        T z=x/a.x,ax2=sqr(a.x);
        auto q=dx-z*a.dx;
        return Make_Hess(z,q/a.x,ddx/a.x-z*a.ddx/a.x-Symmetric_Outer_Product(q,a.dx/ax2));
    }

    AUTO_HESS_EXT operator+(T a) const
    {return Make_Hess(x+a,dx,ddx);}

    AUTO_HESS_EXT operator-(T a) const
    {return Make_Hess(x-a,dx,ddx);}

    auto operator*(T a) const -> decltype(Make_Hess(this->x/a,this->dx/a,this->ddx/a))
    {return Make_Hess(x*a,a*dx,a*ddx);}

    auto operator/(T a) const -> decltype(Make_Hess(this->x/a,this->dx/a,this->ddx/a))
    {return Make_Hess(x/a,dx/a,ddx/a);}

    template<class VEC1,class MAT1>
    void Fill_From(const AUTO_HESS_EXT<T,VEC1,MAT1>& z)
    {x=z.x;::PhysBAM::HETERO_DIFF::Fill_From(dx,z.dx);::PhysBAM::HETERO_DIFF::Fill_From(ddx,z.ddx);}
};

template<class T,class VEC,class MAT> AUTO_HESS_EXT<T,VEC,MAT>
Make_Hess(T x,const GRADIENT<T,VEC>& dx,const HESSIAN<T,MAT>& ddx)
{return AUTO_HESS_EXT<T,VEC,MAT>(x,dx,ddx);}

template<class LAYOUT,class T>
AUTO_HESS_EXT<T,typename EMPTY_VEC<LAYOUT,-1>::TYPE,typename EMPTY_MAT<LAYOUT,-1>::TYPE>
Hess_From_Const(T a)
{return AUTO_HESS_EXT<T,typename EMPTY_VEC<LAYOUT,-1>::TYPE,typename EMPTY_MAT<LAYOUT,-1>::TYPE>(a);}

template<class T,class VEC,class MAT>
inline AUTO_HESS_EXT<T,VEC,MAT> operator+(T c,const AUTO_HESS_EXT<T,VEC,MAT>& a)
{return Make_Hess(c+a.x,a.dx,a.ddx);}

template<class T,class VEC,class MAT>
inline AUTO_HESS_EXT<T,VEC,MAT> operator-(T c,const AUTO_HESS_EXT<T,VEC,MAT>& a)
{return Make_Hess(c-a.x,-a.dx,-a.ddx);}

template<class T,class VEC,class MAT> inline auto
operator*(T c,const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(c*a.x,c*a.dx,c*a.ddx))
{return Make_Hess(c*a.x,c*a.dx,c*a.ddx);}

template<class T,class VEC,class MAT> inline auto
operator/(T c,const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(T(),-T()*a.dx,2*T()/a.x*Outer_Product(a.dx)-T()*a.ddx))
{T z=c/a.x,w=z/a.x;return Make_Hess(z,-w*a.dx,2*w/a.x*Outer_Product(a.dx)-w*a.ddx);}

template<class T,class VEC,class MAT> inline auto
sqrt(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(T(),(a.dx/(2*T())),a.ddx/(2*T())-Outer_Product(a.dx/(2*T()))/T()))
{T s=sqrt(a.x);auto t=a.dx/(2*s);return Make_Hess(s,t,a.ddx/(2*s)-Outer_Product(t)/s);}

template<class T,class VEC,class MAT> inline auto
sqr(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+Outer_Product(a.dx)*2))
{return Make_Hess(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+Outer_Product(a.dx)*2);}

template<class T,class VEC,class MAT> inline auto
cube(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(a.x*T(),3*T()*a.dx,3*T()*a.ddx+Outer_Product(a.dx)*(6*a.x)))
{T sq=sqr(a.x);return Make_Hess(a.x*sq,3*sq*a.dx,3*sq*a.ddx+Outer_Product(a.dx)*(6*a.x));}

template<class T,class VEC,class MAT,class VEC1,class MAT1> inline
AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))>
max(const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC1,MAT1>& b)
{
    AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))> r;
    if(a.x>b.x) r.Fill_From(a);
    else r.Fill_From(b);
    return r;
}

template<class T,class VEC,class MAT> inline auto
max(const AUTO_HESS_EXT<T,VEC,MAT>& a,T b)
    -> decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T())
{
    decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T()) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class T,class VEC,class MAT> inline auto
max(T b,const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T())
{
    decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T()) r;
    if(a.x>b) r.Fill_From(a);
    else r.x=b;
    return r;
}

template<class T,class VEC,class MAT>
inline AUTO_HESS_EXT<T,VEC,MAT>
max(const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC,MAT>& b)
{return a.x>b.x?a:b;}

template<class T,class VEC,class MAT,class VEC1,class MAT1> inline AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))>
min(const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC1,MAT1>& b)
{
    AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1()))> r;
    if(a.x>b.x) r.Fill_From(b);
    else r.Fill_From(a);
    return r;
}

template<class T,class VEC,class MAT> inline auto
min(const AUTO_HESS_EXT<T,VEC,MAT>& a,T b)
    -> decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T())
{
    decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T()) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class T,class VEC,class MAT> inline auto
min(T b,const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T())
{
    decltype(AUTO_HESS_EXT<T,VEC,MAT>()+T()) r;
    if(a.x>b) r.x=b;
    else r.Fill_From(a);
    return r;
}

template<class T,class VEC,class MAT>
inline AUTO_HESS_EXT<T,VEC,MAT>
min(const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC,MAT>& b)
{return a.x>b.x?b:a;}

template<class T,class VEC,class MAT> inline auto
log(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(::std::log(a.x),a.dx/a.x,a.ddx/a.x-Outer_Product(a.dx/a.x)))
{auto z=a.dx/a.x;return Make_Hess(::std::log(a.x),z,a.ddx/a.x-Outer_Product(z));}

template<class T,class VEC,class MAT> inline auto
exp(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*(a.ddx+Outer_Product(a.dx))))
{T s=exp(a.x);return Make_Hess(s,s*a.dx,s*(a.ddx+Outer_Product(a.dx)));}

template<class T,class VEC,class MAT> inline auto
sin(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx-T()*Outer_Product(a.dx)))
{T s=sin(a.x),c=cos(a.x);return Make_Hess(s,c*a.dx,c*a.ddx-s*Outer_Product(a.dx));}

template<class T,class VEC,class MAT> inline auto
cos(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(T(),-T()*a.dx,-T()*a.ddx-T()*Outer_Product(a.dx)))
{T s=sin(a.x),c=cos(a.x);return Make_Hess(c,-s*a.dx,-s*a.ddx-c*Outer_Product(a.dx));}

template<class T,class VEC,class MAT> inline auto
tan(const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx+2*T()*T()*Outer_Product(a.dx)))
{T t=tan(a.x),s=1+t*t;return Make_Hess(t,s*a.dx,s*a.ddx+2*s*t*Outer_Product(a.dx));}

template<class T,class VEC,class MAT> inline AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC()))),decltype(MAT_CHOOSE::Type(MAT(),MAT_NEG::Type(MAT())))>
abs(const AUTO_HESS_EXT<T,VEC,MAT>& a)
{
    AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC()))),decltype(MAT_CHOOSE::Type(MAT(),MAT_NEG::Type(MAT())))> r;
    if(a.x>=0) r.Fill_From(a);
    else r.Fill_From(-a);
    return r;
}

template<class T,class VEC,class MAT,class VEC1,class MAT1> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC1,MAT1>& b)
    -> decltype(Make_Hess(T(),T()*a.dx+T()*b.dx,T()*a.ddx+T()*b.ddx+Outer_Product(T()*b.dx-T()*a.dx)/T()))
{
    T c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return Make_Hess(c,d*a.dx+e*b.dx,d*a.ddx+e*b.ddx+Outer_Product(d*b.dx-e*a.dx)/c);
}

template<class T,class VEC,class MAT> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT>& a,T b)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx+Outer_Product(T()*a.dx)/T()))
{
    T c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return Make_Hess(c,d*a.dx,d*a.ddx+Outer_Product(e*a.dx)/c);
}

template<class T,class VEC,class MAT> inline auto
hypot(T b,const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(hypot(a,b))
{return hypot(a,b);}

template<class T,class VEC,class MAT,class VEC1,class MAT1,class VEC2,class MAT2> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC1,MAT1>& b,const AUTO_HESS_EXT<T,VEC2,MAT2>& c)
    -> decltype(Make_Hess(T(),T()*a.dx+T()*b.dx+T()*c.dx,T()*a.ddx+T()*b.ddx+T()*c.ddx+(Outer_Product(T()*b.dx-T()*a.dx)+Outer_Product(T()*c.dx-T()*b.dx)+Outer_Product(T()*a.dx-T()*c.dx))/T()))
{
    T s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(bb*c.dx-cc*b.dx);
    auto ca=Outer_Product(cc*a.dx-aa*c.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx+cc*c.dx,aa*a.ddx+bb*b.ddx+cc*c.ddx+(ab+bc+ca)/s);
}

template<class T,class VEC,class MAT,class VEC1,class MAT1> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC1,MAT1>& b,T c)
    -> decltype(Make_Hess(T(),T()*a.dx+T()*b.dx,T()*a.ddx+T()*b.ddx+Outer_Product(T()*b.dx-T()*a.dx)/T()+(Outer_Product(b.dx)+Outer_Product(a.dx))*(sqr(T())/T())))
{
    T s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(b.dx);
    auto ca=Outer_Product(a.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx,aa*a.ddx+bb*b.ddx+ab/s+(bc+ca)*(sqr(cc)/s));
}

template<class T,class VEC,class MAT,class VEC1,class MAT1> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT>& a,T c,const AUTO_HESS_EXT<T,VEC1,MAT1>& b)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT,class VEC1,class MAT1> inline auto
hypot(T c,const AUTO_HESS_EXT<T,VEC,MAT>& a,const AUTO_HESS_EXT<T,VEC1,MAT1>& b)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT>& a,T b,T c)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx+Outer_Product(a.dx)*T()))
{
    T t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return Make_Hess(s,d*a.dx,d*a.ddx+Outer_Product(a.dx)*e);
}

template<class T,class VEC,class MAT> inline auto
hypot(T b,const AUTO_HESS_EXT<T,VEC,MAT>& a,T c)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT> inline auto
hypot(T b,T c,const AUTO_HESS_EXT<T,VEC,MAT>& a)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT,class VEC1,class MAT1> inline auto
atan2(const AUTO_HESS_EXT<T,VEC,MAT>& y,const AUTO_HESS_EXT<T,VEC1,MAT1>& x)
    -> decltype(Make_Hess(::std::atan2(y.x,x.x),T()*y.dx-T()*x.dx,T()*y.ddx-T()*x.ddx-Symmetric_Outer_Product(T()*y.dx-T()*x.dx,T()*x.dx+T()*y.dx)))
{
    T c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    auto f=d*y.dx-e*x.dx;
    return Make_Hess(::std::atan2(y.x,x.x),f,d*y.ddx-e*x.ddx-Symmetric_Outer_Product(f,d*x.dx+e*y.dx));
}

template<class T,class VEC,class MAT> inline auto
atan2(const AUTO_HESS_EXT<T,VEC,MAT>& y,T x)
    -> decltype(Make_Hess(::std::atan2(y.x,x),T()*y.dx,T()*y.ddx-Symmetric_Outer_Product(T()*y.dx,T()*y.dx)))
{
    T c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    auto f=d*y.dx;
    return Make_Hess(::std::atan2(y.x,x),f,d*y.ddx-Symmetric_Outer_Product(f,e*y.dx));
}

template<class T,class VEC,class MAT> inline auto
atan2(T y,const AUTO_HESS_EXT<T,VEC,MAT>& x)
    -> decltype(Make_Hess(::std::atan2(y,x.x),-T()*x.dx,-T()*x.ddx-Symmetric_Outer_Product(-T()*x.dx,T()*x.dx)))
{
    T c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    auto f=-e*x.dx;
    return Make_Hess(::std::atan2(y,x.x),f,-e*x.ddx-Symmetric_Outer_Product(f,d*x.dx));
}
template<class TV,class VEC,class MAT> struct AUTO_HESS_EXT_VEC;

template<class TV,class VEC,class MAT> AUTO_HESS_EXT_VEC<TV,VEC,MAT>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const HESSIAN_VEC<TV,MAT>& ddx);

template<class T,class TV,class VEC1,class MAT1,class VEC2,class MAT2> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2>& b)
    -> decltype(Make_Hess_Vec(TV().Cross(TV()),MATRIX<T,TV::m>()*GRADIENT_VEC<TV,VEC2>()-MATRIX<T,TV::m>()*GRADIENT_VEC<TV,VEC1>(),Contract_00(HESSIAN_VEC<TV,MAT1>(),MATRIX<T,TV::m>())-Contract_00(HESSIAN_VEC<TV,MAT2>(),MATRIX<T,TV::m>())+Symmetric_Double_Contract_12_With_Tensor(PERMUTATION_TENSOR<T>(1),GRADIENT_VEC<TV,VEC1>(),GRADIENT_VEC<TV,VEC2>())));

template<class T,class VEC1,class MAT1> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1>& a,const VECTOR<T,3>& b)
    -> decltype(Make_Hess_Vec(VECTOR<T,3>(),-MATRIX<T,3>()*GRADIENT_VEC<VECTOR<T,3>,VEC1>(),Contract_00(HESSIAN_VEC<VECTOR<T,3>,MAT1>(),MATRIX<T,3>())));

template<class TV,class VEC1,class MAT1,class TYPE>
typename enable_if<TV::m!=3,AUTO_HESS_EXT_VEC<TV,VEC1,MAT1> >::type
Cross_Helper(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a,TYPE);

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

    template<class VEC1,class MAT1> auto
    operator+(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const
        -> decltype(Make_Hess_Vec(x+a.x,dx+a.dx,ddx+a.ddx))
    {return Make_Hess_Vec(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class VEC1,class MAT1> auto
    operator-(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const
        -> decltype(Make_Hess_Vec(x-a.x,dx-a.dx,ddx-a.ddx))
    {return Make_Hess_Vec(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class VEC1,class MAT1> auto
    operator*(const AUTO_HESS_EXT<T,VEC1,MAT1>& a) const
        -> decltype(Make_Hess_Vec(this->x*a.x,a.x*this->dx+Outer_Product(this->x,a.dx),this->ddx*a.x+Tensor_Product_0(a.ddx,this->x)+Symmetric_Tensor_Product_12(this->dx,a.dx)))
    {return Make_Hess_Vec(x*a.x,a.x*dx+Outer_Product(x,a.dx),ddx*a.x+Tensor_Product_0(a.ddx,x)+Symmetric_Tensor_Product_12(dx,a.dx));}

    template<class VEC1,class MAT1> auto
    operator/(const AUTO_HESS_EXT<T,VEC1,MAT1>& a) const
        -> decltype(Make_Hess_Vec(this->x/a.x,this->dx/a.x+Outer_Product(x,-a.dx/(a.x*a.x)),this->ddx/a.x+Tensor_Product_0(2/a.x*Outer_Product(a.dx)-a.ddx,this->x/(a.x*a.x))+Symmetric_Tensor_Product_12(this->dx,-a.dx/(a.x*a.x))))
    {auto p=-a.dx/(a.x*a.x);return Make_Hess_Vec(x/a.x,dx/a.x+Outer_Product(x,p),ddx/a.x+Tensor_Product_0(2/a.x*Outer_Product(a.dx)-a.ddx,x/(a.x*a.x))+Symmetric_Tensor_Product_12(dx,p));}

    AUTO_HESS_EXT_VEC operator+(TV a) const
    {return Make_Hess_Vec(x+a,dx,ddx);}

    AUTO_HESS_EXT_VEC operator-(TV a) const
    {return Make_Hess_Vec(x-a,dx,ddx);}

    auto operator*(T a) const -> decltype(Make_Hess_Vec(this->x*a,this->dx*a,this->ddx*a))
    {return Make_Hess_Vec(x*a,dx*a,ddx*a);}

    auto operator/(T a) const -> decltype(Make_Hess_Vec(this->x/a,this->dx/a,this->ddx/a))
    {return Make_Hess_Vec(x/a,dx/a,ddx/a);}

    auto Dot(TV v) const -> decltype(Make_Hess(this->x.Dot(v),this->dx.Transpose_Times(v),Contract_0(this->ddx,v)))
    {return Make_Hess(x.Dot(v),dx.Transpose_Times(v),Contract_0(ddx,v));}

    template<class VEC1,class MAT1> auto
    Dot(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const
        -> decltype(Make_Hess(this->x.Dot(a.x),this->dx.Transpose_Times(a.x)+a.dx.Transpose_Times(this->x),Contract_0(this->ddx,a.x)+Contract_0(a.ddx,this->x)+Symmetric_Transpose_Times(this->dx,a.dx)))
    {return Make_Hess(x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x),Contract_0(ddx,a.x)+Contract_0(a.ddx,x)+Symmetric_Transpose_Times(dx,a.dx));}

    auto Cross(const TV& a) const -> decltype(Cross_Helper(*this,a))
    {return Cross_Helper(*this,a);}

    template<class VEC1,class MAT1>
    auto Cross(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a) const -> decltype(Cross_Helper(*this,a))
    {return Cross_Helper(*this,a);}

    decltype(Make_Hess(x.Magnitude_Squared(),dx.Transpose_Times(x)*2,(Contract_0(ddx,x)+Transpose_Times_Self(dx))*2)) Magnitude_Squared() const
    {return Make_Hess(x.Magnitude_Squared(),dx.Transpose_Times(x)*2,(Contract_0(ddx,x)+Transpose_Times_Self(dx))*2);}

    auto Magnitude() const
        -> decltype(Make_Hess(T(),dx.Transpose_Times(x)/T(),(Contract_0(ddx,x)+Transpose_Times_Self(dx))/T()-Outer_Product(dx.Transpose_Times(x)/T())/T()))
    {T s=x.Magnitude();auto t=dx.Transpose_Times(x)/s;return Make_Hess(s,t,(Contract_0(ddx,x)+Transpose_Times_Self(dx))/s-Outer_Product(t)/s);}
};


template<class T,class TV,class VEC1,class MAT1,class VEC2,class MAT2> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2>& b)
    -> decltype(Make_Hess_Vec(TV().Cross(TV()),MATRIX<T,TV::m>()*GRADIENT_VEC<TV,VEC2>()-MATRIX<T,TV::m>()*GRADIENT_VEC<TV,VEC1>(),Contract_00(HESSIAN_VEC<TV,MAT1>(),MATRIX<T,TV::m>())-Contract_00(HESSIAN_VEC<TV,MAT2>(),MATRIX<T,TV::m>())+Symmetric_Double_Contract_12_With_Tensor(PERMUTATION_TENSOR<T>(1),GRADIENT_VEC<TV,VEC1>(),GRADIENT_VEC<TV,VEC2>())))
{
    MATRIX<T,TV::m> cp_t=MATRIX<T,TV::m>::Cross_Product_Matrix(a.x),cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(b.x);
    return Make_Hess_Vec(a.x.Cross(b.x),cp_t*b.dx-cp_a*a.dx,Contract_00(a.ddx,cp_a)-Contract_00(b.ddx,cp_t)+Symmetric_Double_Contract_12_With_Tensor(PERMUTATION_TENSOR<T>(1),a.dx,b.dx));
}

template<class T,class VEC1,class MAT1> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1>& a,const VECTOR<T,3>& b)
    -> decltype(Make_Hess_Vec(VECTOR<T,3>(),-MATRIX<T,3>()*GRADIENT_VEC<VECTOR<T,3>,VEC1>(),Contract_00(HESSIAN_VEC<VECTOR<T,3>,MAT1>(),MATRIX<T,3>())))
{
    MATRIX<T,3> cp_a=MATRIX<T,3>::Cross_Product_Matrix(b);
    return Make_Hess_Vec(a.x.Cross(b),-cp_a*a.dx,Contract_00(a.ddx,cp_a));
}

template<class TV,class VEC1,class MAT1,class TYPE>
typename enable_if<TV::m!=3,AUTO_HESS_EXT_VEC<TV,VEC1,MAT1> >::type
Cross_Helper(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1>& a,TYPE)
{
    PHYSBAM_FATAL_ERROR("Cross product not defined except in 3D.");
    return a;
}

template<class TV,class VEC,class MAT> inline AUTO_HESS_EXT_VEC<TV,VEC,MAT>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const HESSIAN_VEC<TV,MAT>& ddx)
{return AUTO_HESS_EXT_VEC<TV,VEC,MAT>(x,dx,ddx);}

template<class LAYOUT,int i,class TV> inline
AUTO_HESS_EXT_VEC<TV,typename ONE_NONZERO_VECTOR<i,LAYOUT,TV::m>::TYPE,typename EMPTY_MAT<LAYOUT,TV::m>::TYPE> Hess_From_Var(const TV& v)
{return AUTO_HESS_EXT_VEC<TV,typename ONE_NONZERO_VECTOR<i,LAYOUT,TV::m>::TYPE,typename EMPTY_MAT<LAYOUT,TV::m>::TYPE>(v);}

template<class T,int d,class VEC,class MAT> inline auto
operator*(const AUTO_HESS_EXT<T,VEC,MAT>& a,const VECTOR<T,d>& v)
    -> decltype(Make_Hess_Vec(v*a.x,Outer_Product(v,a.dx),Tensor_Product_0(a.ddx,v)))
{return Make_Hess_Vec(v*a.x,Outer_Product(v,a.dx),Tensor_Product_0(a.ddx,v));}

template<class T,int d,class VEC,class MAT> inline auto
operator*(const VECTOR<T,d>& v,const AUTO_HESS_EXT<T,VEC,MAT>& a) -> decltype(a*v)
{return a*v;}

template<class T,class TV,class VEC,class MAT,class VEC1,class MAT1> inline auto
operator*(const AUTO_HESS_EXT<T,VEC1,MAT1>& a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& v) -> decltype(v*a)
{return v*a;}

template<class TV,class VEC,class MAT> auto
operator*(typename TV::SCALAR a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& v) -> decltype(v*a)
{return v*a;}

template<class TV,class VEC,class MAT> auto
operator+(TV a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& u) -> decltype(Make_Hess_Vec(u.x+a,u.dx,u.ddx))
{return Make_Hess_Vec(u.x+a,u.dx,u.ddx);}

template<class TV,class VEC,class MAT> auto
operator-(TV a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& u)
    -> decltype(Make_Hess_Vec(a-u.x,-u.dx,-u.ddx))
{return Make_Hess_Vec(a-u.x,-u.dx,-u.ddx);}

template<class T,int d,class VEC1,class MAT1> auto
operator/(const VECTOR<T,d>& u,const AUTO_HESS_EXT<T,VEC1,MAT1>& a) -> decltype(u*((T)1/a))
{return u*((T)1/a);}

template<class TV,class VEC,class MAT> auto
operator*(const MATRIX<typename TV::SCALAR,TV::m>& m,const AUTO_HESS_EXT_VEC<TV,VEC,MAT>& v)
    -> decltype(Make_Hess_Vec(m*v.x,m*v.dx,Contract_00(v.ddx,m.Transposed())))
{return Make_Hess_Vec(m*v.x,m*v.dx,Contract_00(v.ddx,m.Transposed()));}
}
using HETERO_DIFF::AUTO_HESS_EXT;
using HETERO_DIFF::AUTO_HESS_EXT_VEC;
using HETERO_DIFF::Hess_From_Const;
using HETERO_DIFF::Hess_From_Var;
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
