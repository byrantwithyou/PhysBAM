//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_HESS
//##################################################################### 
#ifndef __AUTO_HESS__
#define __AUTO_HESS__

#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <cmath>
using ::std::sin;
using ::std::cos;
using ::std::tan;
using ::std::exp;
using ::std::log;
using ::std::abs;
using ::std::sqrt;
using ::PhysBAM::sqr;
using ::PhysBAM::cube;
using ::std::hypot;
namespace PhysBAM{

template<class T,class TV,int Q=3> struct AUTO_HESS;

struct AUTO_HESS_UNUSED
{
    AUTO_HESS_UNUSED x(int)const{return AUTO_HESS_UNUSED();}
    AUTO_HESS_UNUSED operator- () const {return AUTO_HESS_UNUSED();}
    template<class A> AUTO_HESS_UNUSED Transpose_Times(const A&) const{return AUTO_HESS_UNUSED();}
    AUTO_HESS_UNUSED(){}
    template<class T> explicit AUTO_HESS_UNUSED(T){}
    template<class T> AUTO_HESS_UNUSED(T,T){}
    template<class T> AUTO_HESS_UNUSED(T,T,T){}
    AUTO_HESS_UNUSED Transposed() const{return AUTO_HESS_UNUSED();}
    template<class T> void Set_Row(int,T) const{}
    AUTO_HESS_UNUSED Twice_Symmetric_Part() const{return AUTO_HESS_UNUSED();}
    AUTO_HESS_UNUSED Normal_Equations_Matrix() const{return AUTO_HESS_UNUSED();}
    AUTO_HESS_UNUSED Row(int) const{return AUTO_HESS_UNUSED();}
};
template<class T> AUTO_HESS_UNUSED operator+(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED operator-(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED operator*(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED operator/(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED operator+(T,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED operator-(T,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED operator*(T,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
inline AUTO_HESS_UNUSED operator+(AUTO_HESS_UNUSED,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
inline AUTO_HESS_UNUSED operator-(AUTO_HESS_UNUSED,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
inline AUTO_HESS_UNUSED operator*(AUTO_HESS_UNUSED,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<int d,class T> AUTO_HESS_UNUSED Contract(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<int d,int e,class T> AUTO_HESS_UNUSED Contract(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<int d,class T> AUTO_HESS_UNUSED Tensor_Product(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<int d,class T> AUTO_HESS_UNUSED Tensor_Product(T,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<int d> AUTO_HESS_UNUSED Tensor_Product(AUTO_HESS_UNUSED,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED Symmetric_Double_Contract_12(T,AUTO_HESS_UNUSED,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
inline AUTO_HESS_UNUSED Outer_Product(AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED Outer_Product(T,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<int r,int s,class T> AUTO_HESS_UNUSED Symmetric_Tensor_Product(T,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<int r,int s,class T> AUTO_HESS_UNUSED Symmetric_Tensor_Product(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
template<int r,int s> AUTO_HESS_UNUSED Symmetric_Tensor_Product(AUTO_HESS_UNUSED,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED Transpose_Times(T,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<class T> AUTO_HESS_UNUSED Transpose_Times(AUTO_HESS_UNUSED,T){return AUTO_HESS_UNUSED();}
inline AUTO_HESS_UNUSED Transpose_Times(AUTO_HESS_UNUSED,AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<int m,int n> AUTO_HESS_UNUSED Contract(AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<int r,int s> AUTO_HESS_UNUSED Transposed(AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}
template<int r,int s> AUTO_HESS_UNUSED Twice_Symmetric_Part(AUTO_HESS_UNUSED){return AUTO_HESS_UNUSED();}

template<class T,int m,int Q>
struct AUTO_HESS<T,VECTOR<T,m>,Q>
{
    typedef VECTOR<T,m> TV;
    typedef typename std::conditional<(Q&1)!=0,TV,AUTO_HESS_UNUSED>::type DX;
    typedef typename std::conditional<(Q&2)!=0,SYMMETRIC_MATRIX<T,m>,AUTO_HESS_UNUSED>::type DDX;
    
    T x;
    DX dx;
    DDX ddx;

    AUTO_HESS():
        x(T())
    {}

    explicit AUTO_HESS(T x):
        x(x)
    {}

    AUTO_HESS(T x,const DX& dx):
        x(x),dx(dx)
    {}

    AUTO_HESS(T x,const DX& dx,const DDX& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    static AUTO_HESS From_Var(TV v,int i)
    {return AUTO_HESS(v(i),TV::Axis_Vector(i));}

    static AUTO_HESS From_Const(T a)
    {return AUTO_HESS(a,TV());}

    AUTO_HESS operator-() const
    {return AUTO_HESS(-x,-dx,-ddx);}

    AUTO_HESS operator+() const
    {return *this;}

    AUTO_HESS operator+(const AUTO_HESS& a) const
    {return AUTO_HESS(x+a.x,dx+a.dx,ddx+a.ddx);}

    AUTO_HESS operator-(const AUTO_HESS& a) const
    {return AUTO_HESS(x-a.x,dx-a.dx,ddx-a.ddx);}

    AUTO_HESS operator*(const AUTO_HESS& a) const
    {return AUTO_HESS(x*a.x,a.x*dx+x*a.dx,a.x*ddx+x*a.ddx+Outer_Product(dx,a.dx).Twice_Symmetric_Part());}

    AUTO_HESS operator/(const AUTO_HESS& a) const
    {T z=x/a.x,ax2=sqr(a.x);auto q=dx-z*a.dx;return AUTO_HESS(z,q/a.x,ddx/a.x-z*a.ddx/a.x-Outer_Product(q,a.dx/ax2).Twice_Symmetric_Part());}

    AUTO_HESS operator+(T a) const
    {return AUTO_HESS(x+a,dx,ddx);}

    AUTO_HESS operator-(T a) const
    {return AUTO_HESS(x-a,dx,ddx);}

    AUTO_HESS operator*(T a) const
    {return AUTO_HESS(x*a,a*dx,a*ddx);}

    AUTO_HESS operator/(T a) const
    {return AUTO_HESS(x/a,dx/a,ddx/a);}
};

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> operator+(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{return AUTO_HESS<typename TV::SCALAR,TV,Q>(c+a.x,a.dx,a.ddx);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> operator-(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{return AUTO_HESS<typename TV::SCALAR,TV,Q>(c-a.x,-a.dx,-a.ddx);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> operator*(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{return AUTO_HESS<typename TV::SCALAR,TV,Q>(c*a.x,c*a.dx,c*a.ddx);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> operator/(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{typename TV::SCALAR z=c/a.x,w=z/a.x;return AUTO_HESS<typename TV::SCALAR,TV,Q>(z,-w*a.dx,2*w/a.x*Outer_Product(a.dx)-w*a.ddx);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> sqrt(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{typename TV::SCALAR s=sqrt(a.x);auto t=a.dx/(2*s);return AUTO_HESS<typename TV::SCALAR,TV,Q>(s,t,a.ddx/(2*s)-Outer_Product(t)/s);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> sqr(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{return AUTO_HESS<typename TV::SCALAR,TV,Q>(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+Outer_Product(a.dx)*2);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> max(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b)
{return a.x>b.x?a:b;}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> max(typename TV::SCALAR a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b)
{return max(AUTO_HESS<typename TV::SCALAR,TV,Q>(a),b);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> max(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,typename TV::SCALAR b)
{return max(a,AUTO_HESS<typename TV::SCALAR,TV,Q>(b));}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> min(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b)
{return a.x>b.x?b:a;}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> min(typename TV::SCALAR a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b)
{return min(AUTO_HESS<typename TV::SCALAR,TV,Q>(a),b);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> min(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,typename TV::SCALAR b)
{return min(a,AUTO_HESS<typename TV::SCALAR,TV,Q>(b));}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> log(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{auto z=a.dx/a.x;return AUTO_HESS<typename TV::SCALAR,TV,Q>(::std::log(a.x),z,a.ddx/a.x-Outer_Product(z));}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> exp(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{typename TV::SCALAR s=exp(a.x);return AUTO_HESS<typename TV::SCALAR,TV,Q>(s,s*a.dx,s*(a.ddx+Outer_Product(a.dx)));}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> sin(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return AUTO_HESS<typename TV::SCALAR,TV,Q>(s,c*a.dx,c*a.ddx-s*Outer_Product(a.dx));}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> cos(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return AUTO_HESS<typename TV::SCALAR,TV,Q>(c,-s*a.dx,-s*a.ddx-c*Outer_Product(a.dx));}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> tan(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{typename TV::SCALAR t=::std::tan(a.x),s=1+t*t;return AUTO_HESS<typename TV::SCALAR,TV,Q>(t,s*a.dx,s*a.ddx+2*s*t*Outer_Product(a.dx));}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> abs(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{return a.x>=0?a:-a;}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(c,d*a.dx+e*b.dx,d*a.ddx+e*b.ddx+Outer_Product(d*b.dx-e*a.dx)/c);
}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,typename TV::SCALAR b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(c,d*a.dx,d*a.ddx+Outer_Product(e*a.dx)/c);
}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(typename TV::SCALAR b,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{return hypot(a,b);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b,const AUTO_HESS<typename TV::SCALAR,TV,Q>& c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(bb*c.dx-cc*b.dx);
    auto ca=Outer_Product(cc*a.dx-aa*c.dx);
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(s,aa*a.dx+bb*b.dx+cc*c.dx,aa*a.ddx+bb*b.ddx+cc*c.ddx+(ab+bc+ca)/s);
}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b,typename TV::SCALAR c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    auto ab=Outer_Product(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product(b.dx);
    auto ca=Outer_Product(a.dx);
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(s,aa*a.dx+bb*b.dx,aa*a.ddx+bb*b.ddx+ab/s+(bc+ca)*(sqr(cc)/s));
}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b)
{return hypot(a,b,c);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,const AUTO_HESS<typename TV::SCALAR,TV,Q>& b)
{return hypot(a,b,c);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,typename TV::SCALAR b,typename TV::SCALAR c)
{
    typename TV::SCALAR t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(s,d*a.dx,d*a.ddx+Outer_Product(a.dx)*e);
}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(typename TV::SCALAR b,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a,typename TV::SCALAR c)
{return hypot(a,b,c);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> hypot(typename TV::SCALAR b,typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV,Q>& a)
{return hypot(a,b,c);}

inline double hypot(double a,double b,double c)
{return sqrt(a*a+b*b+c*c);}

inline float hypot(float a,float b,float c)
{return sqrt(a*a+b*b+c*c);}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> atan2(const AUTO_HESS<typename TV::SCALAR,TV,Q>& y,const AUTO_HESS<typename TV::SCALAR,TV,Q>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    auto f=d*y.dx-e*x.dx;
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(::std::atan2(y.x,x.x),f,d*y.ddx-e*x.ddx-Outer_Product(f,d*x.dx+e*y.dx).Twice_Symmetric_Part());
}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> atan2(const AUTO_HESS<typename TV::SCALAR,TV,Q>& y,typename TV::SCALAR x)
{
    typename TV::SCALAR c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    auto f=d*y.dx;
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(::std::atan2(y.x,x),f,d*y.ddx-Outer_Product(f,e*y.dx).Twice_Symmetric_Part());
}

template<class TV,int Q>
inline AUTO_HESS<typename TV::SCALAR,TV,Q> atan2(typename TV::SCALAR y,const AUTO_HESS<typename TV::SCALAR,TV,Q>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    auto f=-e*x.dx;
    return AUTO_HESS<typename TV::SCALAR,TV,Q>(::std::atan2(y,x.x),f,-e*x.ddx-Outer_Product(f,d*x.dx).Twice_Symmetric_Part());
}

template<class T,int m,int Q>
AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>
Cross_Product(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& b);

template<class T,int m,int Q>
AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>
Cross_Product(const VECTOR<T,3>& a,const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& b);

template<class T,int m,int Q>
AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>
Cross_Product(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& a,const VECTOR<T,3>& b);

template<class T,int d0,int d1,int m,int Q>
AUTO_HESS<decltype(VECTOR<T,d0>().Cross(VECTOR<T,d1>())),VECTOR<T,m>,Q>
Cross_Product(const AUTO_HESS<VECTOR<T,d0>,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,d1>,VECTOR<T,m>,Q>& b){PHYSBAM_FATAL_ERROR();}

template<class T,int d0,int d1,int m,int Q>
AUTO_HESS<decltype(VECTOR<T,d0>().Cross(VECTOR<T,d1>())),VECTOR<T,m>,Q>
Cross_Product(const VECTOR<T,d0>& a,const AUTO_HESS<VECTOR<T,d1>,VECTOR<T,m>,Q>& b){PHYSBAM_FATAL_ERROR();}

template<class T,int d0,int d1,int m,int Q>
AUTO_HESS<decltype(VECTOR<T,d0>().Cross(VECTOR<T,d1>())),VECTOR<T,m>,Q>
Cross_Product(const AUTO_HESS<VECTOR<T,d0>,VECTOR<T,m>,Q>& a,const VECTOR<T,d1>& b){PHYSBAM_FATAL_ERROR();}

template<class T,int d,int m,int Q>
struct AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
{
    typedef VECTOR<T,m> TV;
    typedef VECTOR<T,d> VEC;
    typedef MATRIX<T,d,TV::m> MAT;
    typedef SYMMETRIC_TENSOR<T,0,d,TV::m> TEN;
    typedef typename std::conditional<(Q&1)!=0,MAT,AUTO_HESS_UNUSED>::type DX;
    typedef typename std::conditional<(Q&2)!=0,TEN,AUTO_HESS_UNUSED>::type DDX;

    VEC x;
    DX dx;
    DDX ddx;

    AUTO_HESS():
        x(VEC())
    {}

    explicit AUTO_HESS(VEC x):
        x(x)
    {}

    AUTO_HESS(VEC x,const DX& dx):
        x(x),dx(dx)
    {}

    AUTO_HESS(VEC x,const DX& dx,const DDX& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    AUTO_HESS(const AUTO_HESS<T,TV,Q>& a)
        :x(a.x),dx(a.dx),ddx(a.ddx)
    {}

    AUTO_HESS(const AUTO_HESS<T,TV,Q>& a,const AUTO_HESS<T,TV,Q>& b)
        :x(a.x,b.x)
    {
        dx.Set_Row(0,a.dx);
        dx.Set_Row(1,b.dx);
        ddx.x(0)=a.ddx;
        ddx.x(1)=b.ddx;
    }

    AUTO_HESS(const AUTO_HESS<T,TV,Q>& a,const AUTO_HESS<T,TV,Q>& b,const AUTO_HESS<T,TV,Q>& c)
        :x(a.x,b.x,c.x)
    {
        ddx.x(0)=a.ddx;
        ddx.x(1)=b.ddx;
        ddx.x(2)=c.ddx;
        dx.Set_Row(0,a.dx);
        dx.Set_Row(1,b.dx);
        dx.Set_Row(2,c.dx);
    }

    static AUTO_HESS From_Var(VEC v)
    {return {v,DX()+1};}

    AUTO_HESS<T,TV,Q> operator()(int i) const
    {return AUTO_HESS<T,TV,Q>(x(i),dx.Row(i),ddx.x(i));}

    void Set_Entry(int i,const AUTO_HESS<T,TV,Q>& a)
    {x(i)=a.x;dx.Set_Row(i,a.dx);ddx.x(i)=a.ddx;}

    AUTO_HESS operator-() const
    {return {-x,-dx,-ddx};}

    AUTO_HESS operator+() const
    {return *this;}

    AUTO_HESS operator+(const AUTO_HESS& a) const
    {return {x+a.x,dx+a.dx,ddx+a.ddx};}

    AUTO_HESS operator-(const AUTO_HESS& a) const
    {return {x-a.x,dx-a.dx,ddx-a.ddx};}

    AUTO_HESS operator*(const AUTO_HESS<T,TV,Q>& a) const
    {return {x*a.x,a.x*dx+Outer_Product(x,a.dx),
        Tensor_Product<0>(a.ddx,x)+Symmetric_Tensor_Product<1,2>(dx,a.dx)+a.x*ddx};}

    AUTO_HESS operator/(const AUTO_HESS<T,TV,Q>& a) const
    {return *this*((T)1/a);}

    AUTO_HESS operator+(VEC a) const
    {return {x+a,dx,ddx};}

    AUTO_HESS operator-(VEC a) const
    {return {x-a,dx,ddx};}

    AUTO_HESS operator*(T a) const
    {return {x*a,a*dx,a*ddx};}

    AUTO_HESS operator/(T a) const
    {return {x/a,dx/a,ddx/a};}

    AUTO_HESS<T,TV,Q> Dot(VEC v) const
    {return {x.Dot(v),dx.Transpose_Times(v),Contract<0>(ddx,v)};}

    AUTO_HESS<T,TV,Q> Dot(const AUTO_HESS& a) const
    {return {x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x),
        Contract<0>(a.ddx,x)+Contract<0>(ddx,a.x)+
        a.dx.Transpose_Times(dx).Twice_Symmetric_Part()};}

    AUTO_HESS<T,TV,Q> Magnitude_Squared() const
    {return Dot(*this);}

    AUTO_HESS<T,TV,Q> Magnitude() const
    {
        T z=x.Magnitude();
        VEC y=x/z;
        auto p=dx.Transpose_Times(y);
        return {z,p,Contract<0>(ddx,y)+(dx.Normal_Equations_Matrix()-Outer_Product(p))/z};
    }

    template<int q>
    auto Cross(const VECTOR<T,q>& a)
    {return Cross_Product(*this,a);}

    template<int q>
    auto Cross(const AUTO_HESS<VECTOR<T,q>,TV,Q>& a)
    {return Cross_Product(*this,a);}
};

template<class T,int m,int Q>
inline AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>
Cross_Product(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& b)
{
    VECTOR<T,3> ab=a.x.Cross(b.x);
    MATRIX<T,3> oa=MATRIX<T,3>::Cross_Product_Matrix(a.x);
    MATRIX<T,3> ob=MATRIX<T,3>::Cross_Product_Matrix(b.x);
    return {ab,oa*b.dx-ob*a.dx,Contract<0,0>(a.ddx,ob)-Contract<0,0>(b.ddx,oa)+Symmetric_Double_Contract_12(PERMUTATION_TENSOR<T>(1),a.dx,b.dx)};
}

template<class T,int m,int Q>
inline AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>
Cross_Product(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& a,const VECTOR<T,3>& b)
{
    VECTOR<T,3> ab=a.x.Cross(b);
    MATRIX<T,3> oa=MATRIX<T,3>::Cross_Product_Matrix(a.x);
    MATRIX<T,3> ob=MATRIX<T,3>::Cross_Product_Matrix(b);
    return {ab,-ob*a.dx,Contract<0,0>(a.ddx,ob)};
}

template<class T,int m,int Q>
inline AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>
Cross_Product(const VECTOR<T,3>& a,const AUTO_HESS<VECTOR<T,3>,VECTOR<T,m>,Q>& b)
{
    VECTOR<T,3> ab=a.Cross(b.x);
    MATRIX<T,3> oa=MATRIX<T,3>::Cross_Product_Matrix(a);
    MATRIX<T,3> ob=MATRIX<T,3>::Cross_Product_Matrix(b.x);
    return {ab,oa*b.dx,-Contract<0,0>(b.ddx,oa)};
}

template<class T,int d,int m,int Q>
inline AUTO_HESS<T,VECTOR<T,m>,Q>
Dot_Product(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{return a.Dot(b);}

template<class T,int d,int m,int Q>
inline AUTO_HESS<T,VECTOR<T,m>,Q>
Dot_Product(const VECTOR<T,d>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{return b.Dot(a);}

template<class T,int d,int m,int Q>
inline AUTO_HESS<T,VECTOR<T,m>,Q>
Dot_Product(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a,const VECTOR<T,d>& b)
{return a.Dot(b);}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator*(const AUTO_HESS<T,VECTOR<T,m>,Q>& a,const VECTOR<T,d>& v)
{return {v*a.x,Outer_Product(v,a.dx),Tensor_Product<0>(a.ddx,v)};}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator+(const VECTOR<T,d>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{return {a+b.x,b.dx,b.ddx};}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator-(const VECTOR<T,d>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{return {a-b.x,-b.dx,-b.ddx};}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator*(const VECTOR<T,d>& v,const AUTO_HESS<T,VECTOR<T,m>,Q>& a)
{return a*v;}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator/(const VECTOR<T,d>& v,const AUTO_HESS<T,VECTOR<T,m>,Q>& a)
{return ((T)1/a)*v;}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator*(const AUTO_HESS<T,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& v)
{return v*a;}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator*(T a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& v)
{return v*a;}

template<class T,int d,int m,int Q>
inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
abs(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{
    AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q> b(a);
    for(int i=0;i<d;i++) if(b.x(i)<0) b.Set_Entry(i,-b(i));
    return b;
}

template<class T,int d,int m,int Q>
inline AUTO_HESS<T,VECTOR<T,m>,Q> max(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{
    AUTO_HESS<T,VECTOR<T,m>,Q> b=a(0);
    for(int i=1;i<d;i++) if(a.x(i)>b.x) b=a(i);
    return b;
}

template<class T,int d,int m,int Q>
inline AUTO_HESS<T,VECTOR<T,m>,Q> min(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{
    AUTO_HESS<T,VECTOR<T,m>,Q> b=a(0);
    for(int i=1;i<d;i++) if(a.x(i)<b.x) b=a(i);
    return b;
}

template<class T> inline VECTOR<T,0> Principal_Curvatures(const AUTO_HESS<T,VECTOR<T,1> >& h)
{return VECTOR<T,0>();}

template<class T> inline VECTOR<T,1> Principal_Curvatures(const AUTO_HESS<T,VECTOR<T,2> >& h)
{VECTOR<T,2> t=h.dx.Perpendicular();return VECTOR<T,1>(t.Dot(h.ddx*t));}

template<class T> VECTOR<T,2> Principal_Curvatures(const AUTO_HESS<T,VECTOR<T,3> >& h);

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
exp(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{VECTOR<T,d> s=exp(a.x);return {s,DIAGONAL_MATRIX<T,d>(s)*a.dx,Contract<0,0>(a.ddx,DIAGONAL_MATRIX<T,d>(s))+Symmetric_Double_Contract_12(DIAGONAL_TENSOR<T,d>(s/2),a.dx,a.dx)};}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
log(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{VECTOR<T,d> b=Inverse(a.x);return {log(a.x),DIAGONAL_MATRIX<T,d>(b)*a.dx,Contract<0,0>(a.ddx,DIAGONAL_MATRIX<T,d>(b))-Symmetric_Double_Contract_12(DIAGONAL_TENSOR<T,d>(0.5*b*b),a.dx,a.dx)};}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
sin(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{VECTOR<T,d> s=sin(a.x),c=cos(a.x);return {s,DIAGONAL_MATRIX<T,d>(c)*a.dx,Contract<0,0>(a.ddx,DIAGONAL_MATRIX<T,d>(c))-Symmetric_Double_Contract_12(DIAGONAL_TENSOR<T,d>(s/2),a.dx,a.dx)};}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
cos(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{VECTOR<T,d> s=sin(a.x),c=cos(a.x);return {c,-DIAGONAL_MATRIX<T,d>(s)*a.dx,Contract<0,0>(a.ddx,-DIAGONAL_MATRIX<T,d>(s))-Symmetric_Double_Contract_12(DIAGONAL_TENSOR<T,d>(c/2),a.dx,a.dx)};}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
Vector_Inverse(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{
    VECTOR<T,d> ia=(T)1/a.x,ia2=ia*ia;
    DIAGONAL_MATRIX<T,d> dia2(-ia2);
    DIAGONAL_TENSOR<T,d> t(ia2*ia);
    return {ia,dia2*a.dx,Contract<0,0>(a.ddx,dia2)+Symmetric_Double_Contract_12(t,a.dx,a.dx)};
}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator*(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{
    DIAGONAL_TENSOR<T,d> o(VECTOR<T,d>()+1);
    DIAGONAL_MATRIX<T,d> oa(a.x),ob(b.x);
    return {a.x*b.x,oa*b.dx+ob*a.dx,Contract<0,0>(a.ddx,ob)+Contract<0,0>(b.ddx,oa)+Symmetric_Double_Contract_12(o,a.dx,b.dx)};
}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator*(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a,const VECTOR<T,d>& b)
{
    DIAGONAL_MATRIX<T,d> ob(b);
    return {a.x*b,ob*a.dx,Contract<0,0>(a.ddx,ob)};
}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator*(const VECTOR<T,d>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{
    DIAGONAL_MATRIX<T,d> oa(a);
    return {a*b.x,oa*b.dx,Contract<0,0>(b.ddx,oa)};
}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator/(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{return a*Vector_Inverse(b);}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator/(const VECTOR<T,d>& a,const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& b)
{return a*Vector_Inverse(b);}

template<class T,int d,int m,int Q> inline AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>
operator/(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a,const VECTOR<T,d>& b)
{return a*((T)1/b);}

template<class T_MAT,class TV,int Q>
struct AUTO_HESS_MAT
{
    STATIC_ASSERT(IS_MATRIX<T_MAT>::value);
    typedef typename TV::SCALAR T;
    typedef T_MAT MAT;
    typedef decltype(Tensor_Product<2>(T_MAT(),TV())) TEN;
    typedef typename std::conditional<Q&1,TEN,AUTO_HESS_UNUSED>::type DX;
    typedef AUTO_HESS_UNUSED DDX;
    enum {m=TV::m,n=T_MAT::m,p=T_MAT::n};

    MAT x;
    DX dx;
    DDX ddx;

    AUTO_HESS_MAT():
        x(MAT())
    {}

    explicit AUTO_HESS_MAT(const MAT& x):
        x(x)
    {}

    AUTO_HESS_MAT(const MAT& x,const DX& dx):
        x(x),dx(dx)
    {}

    AUTO_HESS_MAT(const MAT& x,const DX& dx,const DDX& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    AUTO_HESS<VECTOR<T,p>,VECTOR<T,m>,Q>
    Transpose_Times(AUTO_HESS<VECTOR<T,n>,TV,Q>& a) const
    {
        return {Transpose_Times(x,a.x),Contract<0>(dx,a.x)+Transpose_Times(x,a.dx)};
    }

    AUTO_HESS<VECTOR<T,p>,VECTOR<T,m>,Q>
    Transpose_Times(const VECTOR<T,n>& a) const
    {
        return {::PhysBAM::Transpose_Times(x,a),Contract<0>(dx,a)};
    }

    AUTO_HESS<T,VECTOR<T,m>,Q>
    Trace() const
    {return {x.Trace(),Contract<0,1>(dx)};}

    AUTO_HESS<decltype(x.Transposed()),VECTOR<T,m>,Q>
    Transposed() const
    {return {x.Transposed(),::PhysBAM::Transposed<0,1>(dx)};}

    AUTO_HESS<MATRIX<T,n,p>,VECTOR<T,m>,Q>
    Twice_Symmetric_Part() const
    {return {x.Twice_Symmetric_Part(),::PhysBAM::Twice_Symmetric_Part<0,1>(dx)};}

    template<class T_MAT1>
    typename std::enable_if<IS_MATRIX<T_MAT1>::value,AUTO_HESS<MATRIX<T,p,T_MAT1::n>,VECTOR<T,m>,Q> >::type
    Transpose_Times(const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& a) const
    {return {x.Transpose_Times(a.x),::PhysBAM::Transposed<0,1>(Contract<0,0>(dx,a.x))+Contract<0,0>(a.dx,x)};}

    template<class T_MAT1>
    typename std::enable_if<IS_MATRIX<T_MAT1>::value,AUTO_HESS<MATRIX<T,n,T_MAT1::m>,VECTOR<T,m>,Q> >::type
    Times_Transpose(const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& a) const
    {return {x.Times_Transpose(a.x),Contract<1,1>(dx,a.x)+::PhysBAM::Transposed<0,1>(Contract<1,1>(a.dx,x))};}

    template<class T_MAT1>
    typename std::enable_if<IS_MATRIX<T_MAT1>::value,AUTO_HESS<MATRIX<T,p,T_MAT1::n>,VECTOR<T,m>,Q> >::type
    Transpose_Times(const T_MAT1& a) const
    {return {x.Transpose_Times(a),::PhysBAM::Transposed<0,1>(Contract<0,0>(dx,a))};}

    template<class T_MAT1>
    typename std::enable_if<IS_MATRIX<T_MAT1>::value,AUTO_HESS<MATRIX<T,n,T_MAT1::m>,VECTOR<T,m>,Q> >::type
    Times_Transpose(const T_MAT1& a) const
    {return {x.Times_Transpose(a),Contract<1,1>(dx,a)};}

    AUTO_HESS<SYMMETRIC_MATRIX<T,p>,VECTOR<T,m>,Q>
    Normal_Equations_Matrix() const
    {return {x.Normal_Equations_Matrix(),::PhysBAM::Twice_Symmetric_Part<0,1>(Contract<0,0>(dx,x))};}

    AUTO_HESS<SYMMETRIC_MATRIX<T,n>,VECTOR<T,m>,Q>
    Outer_Product_Matrix() const
    {return {x.Outer_Product_Matrix(),::PhysBAM::Twice_Symmetric_Part<0,1>(Contract<1,1>(dx,x))};}

    template<class T_MAT1>
    typename std::enable_if<IS_MATRIX<T_MAT1>::value,AUTO_HESS<T,VECTOR<T,m>,Q> >::type
    Double_Contract(const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& a) const
    {return {::PhysBAM::Double_Contract(x,a.x),Contract<0,1>(Contract<1,1>(dx,a.x))+Contract<0,1>(Contract<1,1>(a.dx,x))};}

    template<class T_MAT1>
    typename std::enable_if<IS_MATRIX<T_MAT1>::value,AUTO_HESS<T,VECTOR<T,m>,Q> >::type
    Double_Contract(const T_MAT1& a) const
    {return {::PhysBAM::Double_Contract(x,a),Contract<0,1>(Contract<1,1>(dx,a))};}
};

template<class T,int n,int p,int m,int Q>
struct AUTO_HESS<MATRIX<T,n,p>,VECTOR<T,m>,Q>:public AUTO_HESS_MAT<MATRIX<T,n,p>,VECTOR<T,m>,Q>
{
private:
    struct UNUSABLE{};
public:
    typedef AUTO_HESS_MAT<MATRIX<T,n,p>,VECTOR<T,m>,Q> BASE;
    explicit AUTO_HESS(const MATRIX<T,n,p>& x): BASE(x) {}
    AUTO_HESS(const MATRIX<T,n,p>& x,const typename BASE::DX& dx): BASE(x,dx) {}
    AUTO_HESS(const MATRIX<T,n,p>& x,const typename BASE::DX& dx,const typename BASE::DDX& ddx): BASE(x,dx,ddx) {}
    typedef typename std::conditional<(Q&1)&&(n==p),SYMMETRIC_TENSOR<T,2,m,n>,UNUSABLE>::type SYM;
    AUTO_HESS(const MATRIX<T,n,p>& x,const SYM& dx): BASE(x,typename BASE::DX()+dx) {}
    AUTO_HESS(const MATRIX<T,n,p>& x,const SYM& dx,const typename BASE::DDX& ddx): BASE(x,typename BASE::DX()+dx,ddx) {}
};

template<class T,int n,int m,int Q>
struct AUTO_HESS<SYMMETRIC_MATRIX<T,n>,VECTOR<T,m>,Q>:public AUTO_HESS_MAT<SYMMETRIC_MATRIX<T,n>,VECTOR<T,m>,Q>
{
    typedef AUTO_HESS_MAT<SYMMETRIC_MATRIX<T,n>,VECTOR<T,m>,Q> BASE;
    explicit AUTO_HESS(const SYMMETRIC_MATRIX<T,n>& x): BASE(x) {}
    AUTO_HESS(const SYMMETRIC_MATRIX<T,n>& x,const typename BASE::DX& dx): BASE(x,dx) {}
    AUTO_HESS(const SYMMETRIC_MATRIX<T,n>& x,const typename BASE::DX& dx,const typename BASE::DDX& ddx): BASE(x,dx,ddx) {}
};

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a,const AUTO_HESS<T,VECTOR<T,m>,Q>& b)
{
    return {a.x*b.x,b.x*a.dx+Tensor_Product<2>(a.x,b.dx)};
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator*(const T_MAT& a,const AUTO_HESS<T,VECTOR<T,m>,Q>& b)
{
    return {a*b.x,Tensor_Product<2>(a,b.dx)};
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a,const T& b)
{
    return {a.x*b,a.dx*b};
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T,VECTOR<T,m>,Q>& a,const T_MAT& b)
{
    return {b*a.x,Tensor_Product<2>(b,a.dx)};
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T,VECTOR<T,m>,Q>& b,const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a)
{
    return a*b;
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator*(const T& b,const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a)
{
    return a*b;
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator/(const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a,const AUTO_HESS<T,VECTOR<T,m>,Q>& b)
{
    return a*((T)1/b);
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator/(const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a,const T& b)
{
    return a*((T)1/b);
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<T_MAT,VECTOR<T,m>,Q> >::type
operator/(const T_MAT& a,const AUTO_HESS<T,VECTOR<T,m>,Q>& b)
{
    return a*((T)1/b);
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<VECTOR<T,T_MAT::m>,VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a,const AUTO_HESS<VECTOR<T,T_MAT::n>,VECTOR<T,m>,Q>& b)
{
    return {a.x*b.x,Contract<1>(a.dx,b.x)+a.x*b.dx};
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<VECTOR<T,T_MAT::m>,VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T_MAT,VECTOR<T,m>,Q>& a,const VECTOR<T,T_MAT::n>& b)
{
    return {a.x*b,Contract<1>(a.dx,b)};
}

template<class T,class T_MAT,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT>::value,AUTO_HESS<VECTOR<T,T_MAT::m>,VECTOR<T,m>,Q> >::type
operator*(const T_MAT& a,const AUTO_HESS<VECTOR<T,T_MAT::n>,VECTOR<T,m>,Q>& b)
{
    return {a*b.x,a*b.dx,Contract<0,1>(b.ddx,a)};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()+T_MAT1()),VECTOR<T,m>,Q> >::type
operator+(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {a.x+b.x,a.dx+b.dx,a.ddx+b.ddx};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()+T_MAT1()),VECTOR<T,m>,Q> >::type
operator+(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const T_MAT1& b)
{
    return {a.x+b,a.dx,a.ddx};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()+T_MAT1()),VECTOR<T,m>,Q> >::type
operator+(const T_MAT0& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {a+b.x,b.dx,b.ddx};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()-T_MAT1()),VECTOR<T,m>,Q> >::type
operator-(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {a.x-b.x,a.dx-b.dx,a.ddx-b.ddx};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()-T_MAT1()),VECTOR<T,m>,Q> >::type
operator-(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const T_MAT1& b)
{
    return {a.x-b,a.dx,a.ddx};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()-T_MAT1()),VECTOR<T,m>,Q> >::type
operator-(const T_MAT0& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {a-b.x,-b.dx,-b.ddx};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()*T_MAT1()),VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {a.x*b.x,Contract<1,0>(a.dx,b.x)+Contract<0,1>(b.dx,a.x)};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()*T_MAT1()),VECTOR<T,m>,Q> >::type
operator*(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const T_MAT1& b)
{
    return {a.x*b,Contract<1,0>(a.dx,b)};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value,AUTO_HESS<decltype(T_MAT0()*T_MAT1()),VECTOR<T,m>,Q> >::type
operator*(const T_MAT0& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {a*b.x,Contract<0,1>(b.dx,a)};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value&&T_MAT0::m==T_MAT1::m,AUTO_HESS<MATRIX<T,T_MAT0::n,T_MAT1::n>,VECTOR<T,m>,Q> >::type
Transpose_Times(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return a.Transpose_Times(b);
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value&&T_MAT0::m==T_MAT1::m,AUTO_HESS<MATRIX<T,T_MAT0::n,T_MAT1::n>,VECTOR<T,m>,Q> >::type
Transpose_Times(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const T_MAT1& b)
{
    return a.Transpose_Times(b);
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value&&T_MAT0::m==T_MAT1::m,AUTO_HESS<MATRIX<T,T_MAT0::n,T_MAT1::n>,VECTOR<T,m>,Q> >::type
Transpose_Times(const T_MAT0& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {Transpose_Times(a,b.x),Contract<0,0>(b.dx,a)};
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value&&T_MAT0::n==T_MAT1::n,AUTO_HESS<MATRIX<T,T_MAT0::m,T_MAT1::m>,VECTOR<T,m>,Q> >::type
Times_Transpose(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return a.Times_Transpose(b);
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value&&T_MAT0::n==T_MAT1::n,AUTO_HESS<MATRIX<T,T_MAT0::m,T_MAT1::m>,VECTOR<T,m>,Q> >::type
Times_Transpose(const AUTO_HESS<T_MAT0,VECTOR<T,m>,Q>& a,const T_MAT1& b)
{
    return a.Times_Transpose(b);
}

template<class T,class T_MAT0,class T_MAT1,int m,int Q>
typename std::enable_if<IS_MATRIX<T_MAT0>::value&&IS_MATRIX<T_MAT1>::value&&T_MAT0::n==T_MAT1::n,AUTO_HESS<MATRIX<T,T_MAT0::m,T_MAT1::m>,VECTOR<T,m>,Q> >::type
Times_Transpose(const T_MAT0& a,const AUTO_HESS<T_MAT1,VECTOR<T,m>,Q>& b)
{
    return {a*a.x,Contract<1,1>(a.dx,a)};
}

template<class T,int d,int m,int Q>
AUTO_HESS<SYMMETRIC_MATRIX<T,d>,VECTOR<T,m>,Q>
Outer_Product(const AUTO_HESS<VECTOR<T,d>,VECTOR<T,m>,Q>& a)
{
    return {Outer_Product(a.x),Symmetric_Tensor_Product<0,1>(a.dx,a.x)};
}

template<class T,class TV> using AUTO_DIFF=AUTO_HESS<T,TV,1>;
template<class T,class TV> using AUTO_NO_DIFF=AUTO_HESS<T,TV,0>;
}
#endif
