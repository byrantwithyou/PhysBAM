//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_DIFF
//##################################################################### 
#ifndef __AUTO_DIFF__
#define __AUTO_DIFF__

#include <Tools/Matrices/MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

template<class T,class TV> struct AUTO_DIFF;

template<class TV>
struct AUTO_DIFF<typename TV::SCALAR,TV>
{
    typedef typename TV::SCALAR T;

    T x;
    TV dx;
    AUTO_DIFF():
        x(T())
    {}

    AUTO_DIFF(T x,TV dx):
        x(x),dx(dx)
    {}

    static AUTO_DIFF From_Var(TV v,int i)
    {return AUTO_DIFF(v(i),TV::Axis_Vector(i));}

    AUTO_DIFF operator-() const
    {return AUTO_DIFF(-x,-dx);}

    AUTO_DIFF operator+() const
    {return *this;}

    AUTO_DIFF operator+(const AUTO_DIFF& a) const
    {return AUTO_DIFF(x+a.x,dx+a.dx);}

    AUTO_DIFF operator-(const AUTO_DIFF& a) const
    {return AUTO_DIFF(x-a.x,dx-a.dx);}

    AUTO_DIFF operator*(const AUTO_DIFF& a) const
    {return AUTO_DIFF(x*a.x,a.x*dx+x*a.dx);}

    AUTO_DIFF operator/(const AUTO_DIFF& a) const
    {T z=x/a.x,ax2=sqr(a.x);return AUTO_DIFF(z,(dx-z*a.dx)/a.x);}

    AUTO_DIFF operator+(T a) const
    {return AUTO_DIFF(x+a,dx);}

    AUTO_DIFF operator-(T a) const
    {return AUTO_DIFF(x-a,dx);}

    AUTO_DIFF operator*(T a) const
    {return AUTO_DIFF(x*a,a*dx);}

    AUTO_DIFF operator/(T a) const
    {return AUTO_DIFF(x/a,dx/a);}
};

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> operator+(typename TV::SCALAR c,const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_DIFF<typename TV::SCALAR,TV>(c+a.x,a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> operator-(typename TV::SCALAR c,const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_DIFF<typename TV::SCALAR,TV>(c-a.x,-a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> operator*(typename TV::SCALAR c,const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_DIFF<typename TV::SCALAR,TV>(c*a.x,c*a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> operator/(typename TV::SCALAR c,const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR z=c/a.x,w=z/a.x;return AUTO_DIFF<typename TV::SCALAR,TV>(z,-w*a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> sqrt(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=sqrt(a.x);TV t=a.dx/(2*s);return AUTO_DIFF<typename TV::SCALAR,TV>(s,t);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> sqr(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_DIFF<typename TV::SCALAR,TV>(sqr(a.x),2*a.x*a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> max(const AUTO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_DIFF<typename TV::SCALAR,TV>& b)
{return a.x>b.x?a:b;}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> min(const AUTO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_DIFF<typename TV::SCALAR,TV>& b)
{return a.x>b.x?b:a;}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> log(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{TV z=a.dx/a.x;return AUTO_DIFF<typename TV::SCALAR,TV>(::std::log(a.x),z);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> exp(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=exp(a.x);return AUTO_DIFF<typename TV::SCALAR,TV>(s,s*a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> sin(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return AUTO_DIFF<typename TV::SCALAR,TV>(s,c*a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> cos(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return AUTO_DIFF<typename TV::SCALAR,TV>(c,-s*a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> tan(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR t=::std::tan(a.x),s=1+t*t;return AUTO_DIFF<typename TV::SCALAR,TV>(t,s*a.dx);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> abs(const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return a.x>=0?a:-a;}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_DIFF<typename TV::SCALAR,TV>& b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return AUTO_DIFF<typename TV::SCALAR,TV>(c,d*a.dx+e*b.dx);
}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return AUTO_DIFF<typename TV::SCALAR,TV>(c,d*a.dx);
}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return hypot(a,b);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_DIFF<typename TV::SCALAR,TV>& b,const AUTO_DIFF<typename TV::SCALAR,TV>& c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    return AUTO_DIFF<typename TV::SCALAR,TV>(s,aa*a.dx+bb*b.dx+cc*c.dx);
}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_DIFF<typename TV::SCALAR,TV>& b,typename TV::SCALAR c)
{
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    return AUTO_DIFF<typename TV::SCALAR,TV>(s,aa*a.dx+bb*b.dx);
}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR c,const AUTO_DIFF<typename TV::SCALAR,TV>& b)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR c,const AUTO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_DIFF<typename TV::SCALAR,TV>& b)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR b,typename TV::SCALAR c)
{
    typename TV::SCALAR t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return AUTO_DIFF<typename TV::SCALAR,TV>(s,d*a.dx);
}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,const AUTO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR c)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,typename TV::SCALAR c,const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> atan2(const AUTO_DIFF<typename TV::SCALAR,TV>& y,const AUTO_DIFF<typename TV::SCALAR,TV>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    TV f=d*y.dx-e*x.dx;
    return AUTO_DIFF<typename TV::SCALAR,TV>(::std::atan2(y.x,x.x),f);
}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> atan2(const AUTO_DIFF<typename TV::SCALAR,TV>& y,typename TV::SCALAR x)
{
    typename TV::SCALAR c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    TV f=d*y.dx;
    return AUTO_DIFF<typename TV::SCALAR,TV>(::std::atan2(y.x,x),f);
}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> atan2(typename TV::SCALAR y,const AUTO_DIFF<typename TV::SCALAR,TV>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    TV f=-e*x.dx;
    return AUTO_DIFF<typename TV::SCALAR,TV>(::std::atan2(y,x.x),f);
}

template<class TV>
struct AUTO_DIFF<TV,TV>
{
    typedef typename TV::SCALAR T;

    TV x;
    MATRIX<T,TV::m> dx;

    AUTO_DIFF():
        x(T())
    {}

    AUTO_DIFF(TV x,const MATRIX<T,TV::m>& dx):
        x(x),dx(dx)
    {}

    AUTO_DIFF(const AUTO_DIFF<T,TV>& a)
        :x(a.x),dx(a.dx)
    {}

    AUTO_DIFF(const AUTO_DIFF<T,TV>& a,const AUTO_DIFF<T,TV>& b)
        :x(a.x,b.x),dx(a.dx,b.dx)
    {
        dx=dx.Transposed();
    }

    AUTO_DIFF(const AUTO_DIFF<T,TV>& a,const AUTO_DIFF<T,TV>& b,const AUTO_DIFF<T,TV>& c)
        :x(a.x,b.x,c.x),dx(a.dx,b.dx,c.dx)
    {
        dx=dx.Transposed();
    }

    static AUTO_DIFF From_Var(TV v)
    {return AUTO_DIFF(v,MATRIX<T,TV::m>()+1);}

    static AUTO_DIFF From_Const(T a)
    {return AUTO_DIFF(a,TV());}

    AUTO_DIFF<T,TV> operator()(int i) const
    {return AUTO_DIFF<T,TV>(x(i),dx.Row(i));}

    void Set_Entry(int i,const AUTO_DIFF<T,TV>& a)
    {x(i)=a.x;dx.Set_Row(i,a.dx);}

    AUTO_DIFF operator-() const
    {return AUTO_DIFF(-x,-dx);}

    AUTO_DIFF operator+() const
    {return *this;}

    AUTO_DIFF operator+(const AUTO_DIFF& a) const
    {return AUTO_DIFF(x+a.x,dx+a.dx);}

    AUTO_DIFF operator-(const AUTO_DIFF& a) const
    {return AUTO_DIFF(x-a.x,dx-a.dx);}

    AUTO_DIFF operator*(const AUTO_DIFF<T,TV>& a) const
    {return AUTO_DIFF(x*a.x,a.x*dx+MATRIX<T,TV::m>::Outer_Product(x,a.dx));}

    AUTO_DIFF operator/(const AUTO_DIFF<T,TV>& a) const
    {return *this*((T)1/a);}

    AUTO_DIFF operator+(TV a) const
    {return AUTO_DIFF(x+a,dx);}

    AUTO_DIFF operator-(TV a) const
    {return AUTO_DIFF(x-a,dx);}

    AUTO_DIFF operator*(T a) const
    {return AUTO_DIFF(x*a,a*dx);}

    AUTO_DIFF operator/(T a) const
    {return AUTO_DIFF(x/a,dx/a);}

    AUTO_DIFF<T,TV> Dot(TV v) const
    {return AUTO_DIFF<T,TV>(x.Dot(v),dx*v);}

    AUTO_DIFF<T,TV> Dot(const AUTO_DIFF& a) const
    {return AUTO_DIFF<T,TV>(x.Dot(a.x),dx*a.x+a.dx*x);}

    AUTO_DIFF<T,TV> Magnitude_Squared() const
    {return Dot(*this);}

    AUTO_DIFF<T,TV> Magnitude() const;
};

template<class T>
inline AUTO_DIFF<T,VECTOR<T,0> > Magnitude_Helper(const AUTO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_DIFF<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,1> > Magnitude_Helper(const AUTO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return abs(a(0));}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,2> > Magnitude_Helper(const AUTO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return hypot(a(0),a(1));}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,3> > Magnitude_Helper(const AUTO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return hypot(a(0),a(1),a(2));}

template<class TV>
inline AUTO_DIFF<typename TV::SCALAR,TV> AUTO_DIFF<TV,TV>::Magnitude() const {return Magnitude_Helper(*this);}

template<class TV>
inline AUTO_DIFF<TV,TV> operator*(const AUTO_DIFF<typename TV::SCALAR,TV>& a,TV v)
{return AUTO_DIFF<TV,TV>(v*a.x,MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(v,a.dx));}

template<class TV>
inline AUTO_DIFF<TV,TV> operator*(TV v,const AUTO_DIFF<typename TV::SCALAR,TV>& a)
{return a*v;}

template<class TV>
inline AUTO_DIFF<TV,TV> operator*(const AUTO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_DIFF<TV,TV>& v)
{return v*a;}

template<class TV>
inline AUTO_DIFF<TV,TV> operator*(typename TV::SCALAR a,const AUTO_DIFF<TV,TV>& v)
{return v*a;}

template<class T>
inline AUTO_DIFF<VECTOR<T,0>,VECTOR<T,0> > abs(const AUTO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_DIFF<VECTOR<T,0>,VECTOR<T,0> >();}

template<class T>
inline AUTO_DIFF<VECTOR<T,1>,VECTOR<T,1> > abs(const AUTO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return AUTO_DIFF<VECTOR<T,1>,VECTOR<T,1> >(abs(a(0)));}

template<class T>
inline AUTO_DIFF<VECTOR<T,2>,VECTOR<T,2> > abs(const AUTO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return AUTO_DIFF<VECTOR<T,2>,VECTOR<T,2> >(abs(a(0)),abs(a(1)));}

template<class T>
inline AUTO_DIFF<VECTOR<T,3>,VECTOR<T,3> > abs(const AUTO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return AUTO_DIFF<VECTOR<T,3>,VECTOR<T,3> >(abs(a(0)),abs(a(1)),abs(a(2)));}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,0> > max(const AUTO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_DIFF<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,1> > max(const AUTO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,2> > max(const AUTO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return max(a(0),a(1));}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,3> > max(const AUTO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return max(a(0),max(a(1),a(2)));}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,0> > min(const AUTO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_DIFF<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,1> > min(const AUTO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,2> > min(const AUTO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return min(a(0),a(1));}

template<class T>
inline AUTO_DIFF<T,VECTOR<T,3> > min(const AUTO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return min(a(0),min(a(1),a(2)));}

}
#endif
