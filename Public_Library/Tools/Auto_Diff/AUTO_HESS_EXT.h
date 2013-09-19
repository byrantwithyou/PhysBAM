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

template<class TV,class G,class H> struct AUTO_HESS_EXT;

#define MK_AH(TV,GFLAGS,HFLAGS) AUTO_HESS_EXT<TV,decltype(GFLAGS),decltype(HFLAGS)>

template<class TV,class G,class H> AUTO_HESS_EXT<TV,G,H>
Make_Hess(typename TV::SCALAR x,const GRADIENT<TV,G>& dx,const HESSIAN<TV,H>& ddx);

template<class TV,class G,class H>
struct AUTO_HESS_EXT
{
    typedef typename TV::SCALAR T;

    T x;
    GRADIENT<TV,G> dx;
    HESSIAN<TV,H> ddx;

    AUTO_HESS_EXT(T z=T()): x(z)
    {}

    AUTO_HESS_EXT(T z,const GRADIENT<TV,G>& dz,const HESSIAN<TV,H>& ddz)
        :x(z),dx(dz),ddx(ddz)
    {}

    MK_AH(TV,-G(),-H()) operator-() const
    {return Make_Hess(-x,-dx,-ddx);}

    AUTO_HESS_EXT operator+() const
    {return *this;}

    template<class G2,class H2> MK_AH(TV,G()+G2(),H()+H2())
    operator+(const AUTO_HESS_EXT<TV,G2,H2>& a) const
    {return Make_Hess(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class G2,class H2> MK_AH(TV,G()-G2(),H()-H2())
    operator-(const AUTO_HESS_EXT<TV,G2,H2>& a) const
    {return Make_Hess(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class G2,class H2> MK_AH(TV,G()+G2(),H()*1+H2()*1+Sym_Outer(G(),G2()))
    operator*(const AUTO_HESS_EXT<TV,G2,H2>& a) const
    {return Make_Hess(x*a.x,a.x*dx+x*a.dx,a.x*ddx+x*a.ddx+Symmetric_Outer_Product(dx,a.dx));}

    template<class G2,class H2> MK_AH(TV,G()+G2(),H()*1+H2()*1+Sym_Outer(G()-G2(),G2()))
    operator/(const AUTO_HESS_EXT<TV,G2,H2>& a) const
    {T z=x/a.x,ax2=sqr(a.x);auto q=dx-z*a.dx;return Make_Hess(z,q/a.x,ddx/a.x-z*a.ddx/a.x-Symmetric_Outer_Product(q,a.dx/ax2));}

    AUTO_HESS_EXT operator+(T a) const
    {return Make_Hess(x+a,dx,ddx);}

    AUTO_HESS_EXT operator-(T a) const
    {return Make_Hess(x-a,dx,ddx);}

    MK_AH(TV,G()*1,H()*1) operator*(T a) const
    {return Make_Hess(x*a,a*dx,a*ddx);}

    MK_AH(TV,G()*1,H()*1) operator/(T a) const
    {return Make_Hess(x/a,dx/a,ddx/a);}

    template<class G2,class H2>
    void Fill_From(const AUTO_HESS_EXT<TV,G2,H2>& z)
    {x=z.x;dx.Fill_From(z.dx);ddx.Fill_From(z.ddx);}
};

template<class TV,class G,class H> AUTO_HESS_EXT<TV,G,H>
Make_Hess(typename TV::SCALAR x,const GRADIENT<TV,G>& dx,const HESSIAN<TV,H>& ddx)
{return AUTO_HESS_EXT<TV,G,H>(x,dx,ddx);}

template<int n,int i,class TV>
AUTO_HESS_EXT<TV,CV_FLAGS<n,(1<<i)>,MM_FLAGS<n,0,0> > From_Var(TV v,int j)
{AUTO_HESS_EXT<TV,CV_FLAGS<n,(1<<i)>,MM_FLAGS<n,0,0> > r(v(j));r.dx.template Set_Entry<i>(TV::Axis_Vector(j));return r;}

template<class TV,int n>
AUTO_HESS_EXT<TV,CV_FLAGS<n,0>,MM_FLAGS<n,0,0> > From_Const(typename TV::SCALAR a)
{return AUTO_HESS_EXT<TV,CV_FLAGS<n,0>,MM_FLAGS<n,0,0> >(a);}

template<class TV,class G,class H>
inline AUTO_HESS_EXT<TV,G,H> operator+(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,G,H>& a)
{return Make_Hess(c+a.x,a.dx,a.ddx);}

template<class TV,class G,class H>
inline AUTO_HESS_EXT<TV,G,H> operator-(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,G,H>& a)
{return Make_Hess(c-a.x,-a.dx,-a.ddx);}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,H()*1)
operator*(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,G,H>& a)
{return Make_Hess(c*a.x,c*a.dx,c*a.ddx);}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
operator/(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,G,H>& a)
{typename TV::SCALAR z=c/a.x,w=z/a.x;return Make_Hess(z,-w*a.dx,2*w/a.x*Outer_Product(a.dx)-w*a.ddx);}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
sqrt(const AUTO_HESS_EXT<TV,G,H>& a)
{typename TV::SCALAR s=sqrt(a.x);auto t=a.dx/(2*s);return Make_Hess(s,t,a.ddx/(2*s)-Outer_Product(t)/s);}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
sqr(const AUTO_HESS_EXT<TV,G,H>& a)
{return Make_Hess(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+Outer_Product(a.dx)*2);}

template<class TV,class G,class H,class G2,class H2> inline MK_AH(TV,Choose(G(),G2()),Choose(H(),H2()))
max(const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G2,H2>& b)
{
    MK_AH(TV,Choose(G(),G2()),Choose(H(),H2())) r;
    if(a.x>b.x) r.Fill_From(a);
    else r.Fill_From(b);
    return r;
}

template<class TV,class G,class H>
inline AUTO_HESS_EXT<TV,G,H>
max(const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G,H>& b)
{return a.x>b.x?a:b;}

template<class TV,class G,class H,class G2,class H2> inline MK_AH(TV,Choose(G(),G2()),Choose(H(),H2()))
min(const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G2,H2>& b)
{
    MK_AH(TV,Choose(G(),G2()),Choose(H(),H2())) r;
    if(a.x>b.x) r.Fill_From(b);
    else r.Fill_From(a);
    return r;
}

template<class TV,class G,class H>
inline AUTO_HESS_EXT<TV,G,H>
min(const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G,H>& b)
{return a.x>b.x?b:a;}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
log(const AUTO_HESS_EXT<TV,G,H>& a)
{auto z=a.dx/a.x;return Make_Hess(::std::log(a.x),z,a.ddx/a.x-Outer_Product(z));}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
exp(const AUTO_HESS_EXT<TV,G,H>& a)
{typename TV::SCALAR s=exp(a.x);return Make_Hess(s,s*a.dx,s*(a.ddx+Outer_Product(a.dx)));}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
sin(const AUTO_HESS_EXT<TV,G,H>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return Make_Hess(s,c*a.dx,c*a.ddx-s*Outer_Product(a.dx));}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
cos(const AUTO_HESS_EXT<TV,G,H>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return Make_Hess(c,-s*a.dx,-s*a.ddx-c*Outer_Product(a.dx));}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
tan(const AUTO_HESS_EXT<TV,G,H>& a)
{typename TV::SCALAR t=tan(a.x),s=1+t*t;return Make_Hess(t,s*a.dx,s*a.ddx+2*s*t*Outer_Product(a.dx));}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,H()*1)
abs(const AUTO_HESS_EXT<TV,G,H>& a)
{
    MK_AH(TV,G()*1,H()*1) r;
    if(a.x>=0) r.Fill_From(a);
    else r.Fill_From(-a);
    return r;
}



template<class TV,class G,class H,class G2,class H2> MK_AH(TV,G()+G2(),H()*1+H2()*1+Outer(G()+G2()))
hypot(const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G2,H2>& b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return Make_Hess(c,d*a.dx+e*b.dx,d*a.ddx+e*b.ddx+Outer_Product(d*b.dx-e*a.dx)/c);
}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
hypot(const AUTO_HESS_EXT<TV,G,H>& a,typename TV::SCALAR b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return Make_Hess(c,d*a.dx,d*a.ddx+Outer_Product(e*a.dx)/c);
}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
hypot(typename TV::SCALAR b,const AUTO_HESS_EXT<TV,G,H>& a)
{return hypot(a,b);}

template<class TV,class G,class H,class G2,class H2,class G3,class H3>
MK_AH(TV,G()+G2()+G3(),1*H()+1*H2()+1*H3()+Outer(G()+G2())+Outer(G()+G3())+Outer(G2()+G3()))
hypot(const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G2,H2>& b,const AUTO_HESS_EXT<TV,G3,H3>& c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(bb*c.dx-cc*b.dx);
    auto ca=Outer_Product(cc*a.dx-aa*c.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx+cc*c.dx,aa*a.ddx+bb*b.ddx+cc*c.ddx+(ab+bc+ca)/s);
}

template<class TV,class G,class H,class G2,class H2> MK_AH(TV,G()+G2(),H()*1+H2()*1+Outer(G()+G2()))
hypot(const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G2,H2>& b,typename TV::SCALAR c)
{
    typedef SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> SM;
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    SM ab=Outer_Product(aa*b.dx-bb*a.dx);
    SM bc=Outer_Product(b.dx);
    SM ca=Outer_Product(a.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx,aa*a.ddx+bb*b.ddx+ab/s+(bc+ca)*(sqr(cc)/s));
}

template<class TV,class G,class H,class G2,class H2> MK_AH(TV,G()+G2(),H()*1+H2()*1+Outer(G()+G2()))
hypot(const AUTO_HESS_EXT<TV,G,H>& a,typename TV::SCALAR c,const AUTO_HESS_EXT<TV,G2,H2>& b)
{return hypot(a,b,c);}

template<class TV,class G,class H,class G2,class H2> MK_AH(TV,G()+G2(),H()*1+H2()*1+Outer(G()+G2()))
hypot(typename TV::SCALAR c,const AUTO_HESS_EXT<TV,G,H>& a,const AUTO_HESS_EXT<TV,G2,H2>& b)
{return hypot(a,b,c);}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
hypot(const AUTO_HESS_EXT<TV,G,H>& a,typename TV::SCALAR b,typename TV::SCALAR c)
{
    typename TV::SCALAR t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return Make_Hess(s,d*a.dx,d*a.ddx+Outer_Product(a.dx)*e);
}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
hypot(typename TV::SCALAR b,const AUTO_HESS_EXT<TV,G,H>& a,typename TV::SCALAR c)
{return hypot(a,b,c);}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
hypot(typename TV::SCALAR b,typename TV::SCALAR c,const AUTO_HESS_EXT<TV,G,H>& a)
{return hypot(a,b,c);}

template<class TV,class G,class H,class G2,class H2> MK_AH(TV,G()+G2(),H()*1+H2()*1+Outer(G()+G2()))
atan2(const AUTO_HESS_EXT<TV,G,H>& y,const AUTO_HESS_EXT<TV,G2,H2>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    auto f=d*y.dx-e*x.dx;
    return Make_Hess(::std::atan2(y.x,x.x),f,d*y.ddx-e*x.ddx-Symmetric_Outer_Product(f,d*x.dx+e*y.dx));
}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
atan2(const AUTO_HESS_EXT<TV,G,H>& y,typename TV::SCALAR x)
{
    typename TV::SCALAR c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    auto f=d*y.dx;
    return Make_Hess(::std::atan2(y.x,x),f,d*y.ddx-Symmetric_Outer_Product(f,e*y.dx));
}

template<class TV,class G,class H> inline MK_AH(TV,G()*1,Outer(G())+H()*1)
atan2(typename TV::SCALAR y,const AUTO_HESS_EXT<TV,G,H>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    auto f=-e*x.dx;
    return Make_Hess(::std::atan2(y,x.x),f,-e*x.ddx-Symmetric_Outer_Product(f,d*x.dx));
}
template<class TV,class H,class R> struct AUTO_HESS_EXT_VEC;

template<class TV,class H,class R> AUTO_HESS_EXT_VEC<TV,H,R>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,H>& dx,const HESSIAN_VEC<TV,R>& ddx);

#define MK_AHV(TV,HFLAGS,RFLAGS) AUTO_HESS_EXT_VEC<TV,decltype(HFLAGS),decltype(RFLAGS)>

template<class TV,class H,class R,class H2,class R2>
struct DOT_HELPER
{
    typedef MK_AH(TV,H()*V_FLAGS<1>()+H2()*V_FLAGS<1>(),C0(R(),V_FLAGS<1>())+C0(R2(),V_FLAGS<1>())+Stt(H(),H2())) TYPE;
};

template<class TV,class H,class R,class H2,class R2>
struct CROSS_HELPER
{
    typedef MK_AHV(TV,(H2()*M_FLAGS<1,0>()-H()*M_FLAGS<1,0>()),(C0(R(),M_FLAGS<1,0>())+C0(R2(),M_FLAGS<1,0>())+typename SYMMETRIC_DOUBLE_CONTRACT_HELPER<T_FLAGS<0,1,1>,H,H2>::TYPE())) TYPE;
};

template<class TV,class H,class R>
struct AUTO_HESS_EXT_VEC
{
    typedef typename TV::SCALAR T;

    TV x;
    GRADIENT_VEC<TV,H> dx;
    HESSIAN_VEC<TV,R> ddx;

    AUTO_HESS_EXT_VEC(TV x=TV()):
        x(x)
    {}

    AUTO_HESS_EXT_VEC(TV x,const GRADIENT_VEC<TV,H>& dx,const HESSIAN_VEC<TV,R>& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    MK_AHV(TV,-H(),-R()) operator-() const
    {return Make_Hess_Vec(-x,-dx,-ddx);}

    AUTO_HESS_EXT_VEC operator+() const
    {return *this;}

    template<class H2,class R2>
    MK_AHV(TV,H()+H2(),R()+R2()) operator+(const AUTO_HESS_EXT_VEC<TV,H2,R2>& a) const
    {return Make_Hess_Vec(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class H2,class R2>
    MK_AHV(TV,H()-H2(),R()-R2()) operator-(const AUTO_HESS_EXT_VEC<TV,H2,R2>& a) const
    {return Make_Hess_Vec(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class G2,class H2>
    MK_AHV(TV,H()*1+Outer(V_FLAGS<1>(),G2()),R()+Tp0(V_FLAGS<1>(),H2())+Stp12(G2(),H()))
    operator*(const AUTO_HESS_EXT<TV,G2,H2>& a) const
    {return Make_Hess_Vec(x*a.x,a.x*dx+Outer_Product(x,a.dx),ddx*a.x+Tensor_Product_0(a.ddx,x)+Symmetric_Tensor_Product_12(dx,a.dx));}

    template<class G2,class H2>
    MK_AHV(TV,H()*1+Outer(V_FLAGS<1>(),G2()*1),R()+Tp0(V_FLAGS<1>(),Outer(G2())+H2()*1)+Stp12(G2()*1,H()))
    operator/(const AUTO_HESS_EXT<TV,G2,H2>& a) const
    {return *this*((T)1/a);}

    AUTO_HESS_EXT_VEC operator+(TV a) const
    {return Make_Hess_Vec(x+a,dx,ddx);}

    AUTO_HESS_EXT_VEC operator-(TV a) const
    {return Make_Hess_Vec(x-a,dx,ddx);}

    MK_AHV(TV,H()*1,R()*1) operator*(T a) const
    {return Make_Hess_Vec(x*a,dx*a,ddx*a);}

    MK_AHV(TV,H()*1,R()*1) operator/(T a) const
    {return Make_Hess_Vec(x/a,dx/a,ddx/a);}

    MK_AH(TV,H()*V_FLAGS<1>(),C0(R(),V_FLAGS<1>())) Dot(TV v) const
    {return Make_Hess(x.Dot(v),dx.Transpose_Times(v),Contract_0(ddx,v));}

    template<class H2,class R2>
    typename DOT_HELPER<TV,H,R,H2,R2>::TYPE
    Dot(const AUTO_HESS_EXT_VEC<TV,H2,R2>& a) const
    {return Make_Hess(x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x),Contract_0(ddx,a.x)+Contract_0(a.ddx,x)+Symmetric_Transpose_Times(dx,a.dx));}

    template<class H2,class R2>
    typename CROSS_HELPER<TV,H,R,H2,R2>::TYPE
    Cross(const AUTO_HESS_EXT_VEC<TV,H2,R2>& a) const
    {
        MATRIX<T,TV::m> cp_t=MATRIX<T,TV::m>::Cross_Product_Matrix(x),cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(a.x);
        return Make_Hess_Vec(x.Cross(a.x),cp_t*a.dx-cp_a*dx,Contract_0(ddx,cp_a)-Contract_0(a.ddx,cp_t)+Symmetric_Double_Contract_12_With_Tensor(PERM_TENSOR<TV>(1),dx,a.dx));
    }

    MK_AH(TV,H()*V_FLAGS<1>(),C0(R(),V_FLAGS<1>())+Stt(H(),H())) Magnitude_Squared() const
    {return Make_Hess(x.Magnitude_Squared(),dx.Transpose_Times(x)*2,(Contract_0(ddx,x)+Transpose_Times_Self(dx))*2);}

    decltype(sqrt(MK_AH(TV,H()*V_FLAGS<1>(),C0(R(),V_FLAGS<1>())+Stt(H(),H()))())) Magnitude() const{return sqrt(Magnitude_Squared());}
};

template<class TV,class H,class R> inline AUTO_HESS_EXT_VEC<TV,H,R>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,H>& dx,const HESSIAN_VEC<TV,R>& ddx)
{return AUTO_HESS_EXT_VEC<TV,H,R>(x,dx,ddx);}

template<int n,int i,class TV> inline
AUTO_HESS_EXT_VEC<TV,CM_FLAGS<n,0,(1<<i)>,MT_FLAGS<n,0,0,0> > From_Var(TV v)
{return AUTO_HESS_EXT_VEC<TV,CM_FLAGS<n,0,(1<<i)>,MT_FLAGS<n,0,0,0> >(v);}

template<class TV,class G,class H> inline MK_AHV(TV,Outer(V_FLAGS<1>(),G()),Tp0(V_FLAGS<1>(),H()))
operator*(const AUTO_HESS_EXT<TV,G,H>& a,TV v)
{return Make_Hess_Vec(v*a.x,Outer_Product(v,a.dx),Tensor_Product_0(a.ddx,v));}

template<class TV,class G,class H> inline MK_AHV(TV,Outer(V_FLAGS<1>(),G()),Tp0(V_FLAGS<1>(),H()))
operator*(TV v,const AUTO_HESS_EXT<TV,G,H>& a)
{return a*v;}

template<class TV,class H,class R,class G2,class H2> inline
MK_AHV(TV,H()*1+Outer(V_FLAGS<1>(),G2()),R()+Tp0(V_FLAGS<1>(),H2())+Stp12(G2(),H()))
operator*(const AUTO_HESS_EXT<TV,G2,H2>& a,const AUTO_HESS_EXT_VEC<TV,H,R>& v)
{return v*a;}

template<class TV,class H,class R> MK_AHV(TV,H()*1,R()*1)
operator*(typename TV::SCALAR a,const AUTO_HESS_EXT_VEC<TV,H,R>& v)
{return v*a;}

template<class TV,class H,class R> AUTO_HESS_EXT_VEC<TV,H,R>
operator+(TV a,const AUTO_HESS_EXT_VEC<TV,H,R>& u)
{return Make_Hess_Vec(u.x+a,u.dx,u.ddx);}

template<class TV,class H,class R> MK_AHV(TV,-H(),-R())
operator-(TV a,const AUTO_HESS_EXT_VEC<TV,H,R>& u)
{return Make_Hess_Vec(a-u.x,-u.dx,-u.ddx);}

template<class TV,class G2,class H2>
MK_AHV(TV,Outer(V_FLAGS<1>(),G2()*1),Tp0(V_FLAGS<1>(),Outer(G2())+H2()*1))
operator/(TV u,const AUTO_HESS_EXT<TV,G2,H2>& a)
{return u*((typename TV::SCALAR)1/a);}

#if 0
template<class T>
inline AUTO_HESS_EXT<VECTOR<T,0>,VECTOR<T,0> > abs(const AUTO_HESS_EXT<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_HESS_EXT<VECTOR<T,0>,VECTOR<T,0> >();}

template<class T>
inline AUTO_HESS_EXT<VECTOR<T,1>,VECTOR<T,1> > abs(const AUTO_HESS_EXT<VECTOR<T,1>,VECTOR<T,1> >& a)
{return AUTO_HESS_EXT<VECTOR<T,1>,VECTOR<T,1> >(abs(a(0)));}

template<class T>
inline AUTO_HESS_EXT<VECTOR<T,2>,VECTOR<T,2> > abs(const AUTO_HESS_EXT<VECTOR<T,2>,VECTOR<T,2> >& a)
{return AUTO_HESS_EXT<VECTOR<T,2>,VECTOR<T,2> >(abs(a(0)),abs(a(1)));}

template<class T>
inline AUTO_HESS_EXT<VECTOR<T,3>,VECTOR<T,3> > abs(const AUTO_HESS_EXT<VECTOR<T,3>,VECTOR<T,3> >& a)
{return AUTO_HESS_EXT<VECTOR<T,3>,VECTOR<T,3> >(abs(a(0)),abs(a(1)),abs(a(2)));}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,0> > max(const AUTO_HESS_EXT<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_HESS_EXT<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,1> > max(const AUTO_HESS_EXT<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,2> > max(const AUTO_HESS_EXT<VECTOR<T,2>,VECTOR<T,2> >& a)
{return max(a(0),a(1));}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,3> > max(const AUTO_HESS_EXT<VECTOR<T,3>,VECTOR<T,3> >& a)
{return max(a(0),max(a(1),a(2)));}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,0> > min(const AUTO_HESS_EXT<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_HESS_EXT<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,1> > min(const AUTO_HESS_EXT<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,2> > min(const AUTO_HESS_EXT<VECTOR<T,2>,VECTOR<T,2> >& a)
{return min(a(0),a(1));}

template<class T>
inline AUTO_HESS_EXT<T,VECTOR<T,3> > min(const AUTO_HESS_EXT<VECTOR<T,3>,VECTOR<T,3> >& a)
{return min(a(0),min(a(1),a(2)));}
#endif
}
using HETERO_DIFF::AUTO_HESS_EXT;
using HETERO_DIFF::AUTO_HESS_EXT_VEC;
using HETERO_DIFF::From_Const;
using HETERO_DIFF::From_Var;
using HETERO_DIFF::sin;
using HETERO_DIFF::cos;
using HETERO_DIFF::tan;
using HETERO_DIFF::atan2;
using HETERO_DIFF::exp;
using HETERO_DIFF::log;
using HETERO_DIFF::min;
using HETERO_DIFF::max;
using HETERO_DIFF::abs;
using HETERO_DIFF::hypot;
using HETERO_DIFF::sqr;
using HETERO_DIFF::sqrt;
}
#endif
