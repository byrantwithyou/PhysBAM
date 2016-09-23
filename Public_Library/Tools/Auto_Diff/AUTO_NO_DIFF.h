//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_NO_DIFF
//##################################################################### 
#ifndef __AUTO_NO_DIFF__
#define __AUTO_NO_DIFF__

#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

template<class T,class TV> struct AUTO_NO_DIFF;

template<class TV>
struct AUTO_NO_DIFF<typename TV::SCALAR,TV>
{
    typedef typename TV::SCALAR T;

    T x;
    AUTO_NO_DIFF():
        x(T())
    {}

    AUTO_NO_DIFF(T x):
        x(x)
    {}

    static AUTO_NO_DIFF From_Var(TV v,int i)
    {return AUTO_NO_DIFF(v(i),TV::Axis_Vector(i));}

    AUTO_NO_DIFF operator-() const
    {return AUTO_NO_DIFF(-x);}

    AUTO_NO_DIFF operator+() const
    {return *this;}

    AUTO_NO_DIFF operator+(const AUTO_NO_DIFF& a) const
    {return AUTO_NO_DIFF(x+a.x);}

    AUTO_NO_DIFF operator-(const AUTO_NO_DIFF& a) const
    {return AUTO_NO_DIFF(x-a.x);}

    AUTO_NO_DIFF operator*(const AUTO_NO_DIFF& a) const
    {return AUTO_NO_DIFF(x*a.x);}

    AUTO_NO_DIFF operator/(const AUTO_NO_DIFF& a) const
    {return AUTO_NO_DIFF(x/a.x);}

    AUTO_NO_DIFF operator+(T a) const
    {return AUTO_NO_DIFF(x+a);}

    AUTO_NO_DIFF operator-(T a) const
    {return AUTO_NO_DIFF(x-a);}

    AUTO_NO_DIFF operator*(T a) const
    {return AUTO_NO_DIFF(x*a);}

    AUTO_NO_DIFF operator/(T a) const
    {return AUTO_NO_DIFF(x/a);}
};

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> operator+(typename TV::SCALAR c,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(c+a.x);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> operator-(typename TV::SCALAR c,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(c-a.x);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> operator*(typename TV::SCALAR c,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(c*a.x);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> operator/(typename TV::SCALAR c,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(c/a.x);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> sqrt(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sqrt(a.x));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> sqr(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sqr(a.x));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> max(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& b)
{return a.x>b.x?a:b;}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> min(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& b)
{return a.x>b.x?b:a;}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> log(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(::std::log(a.x));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> exp(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(exp(a.x));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> sin(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sin(a.x));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> cos(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(cos(a.x));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> tan(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return AUTO_NO_DIFF<typename TV::SCALAR,TV>(::std::tan(a.x));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> abs(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return a.x>=0?a:-a;}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& b)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sqrt(sqr(a.x)+sqr(b.x)));
}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR b)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sqrt(sqr(a.x)+sqr(b)));
}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return hypot(a,b);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& b,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& c)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)));
}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& b,typename TV::SCALAR c)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sqrt(sqr(a.x)+sqr(b.x)+sqr(c)));
}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR c,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& b)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR c,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& b)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR b,typename TV::SCALAR c)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(sqrt(sqr(a.x)+sqr(b)+sqr(c)));
}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,typename TV::SCALAR c)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,typename TV::SCALAR c,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> atan2(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& y,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& x)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(::std::atan2(y.x,x.x));
}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> atan2(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& y,typename TV::SCALAR x)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(::std::atan2(y.x,x));
}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> atan2(typename TV::SCALAR y,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& x)
{
    return AUTO_NO_DIFF<typename TV::SCALAR,TV>(::std::atan2(y,x.x));
}

template<class TV>
struct AUTO_NO_DIFF<TV,TV>
{
    typedef typename TV::SCALAR T;

    TV x;

    AUTO_NO_DIFF():
        x(T())
    {}

    AUTO_NO_DIFF(TV x):
        x(x)
    {}

    AUTO_NO_DIFF(const AUTO_NO_DIFF<T,TV>& a)
        :x(a.x)
    {}

    AUTO_NO_DIFF(const AUTO_NO_DIFF<T,TV>& a,const AUTO_NO_DIFF<T,TV>& b)
        :x(a.x,b.x)
    {
    }

    AUTO_NO_DIFF(const AUTO_NO_DIFF<T,TV>& a,const AUTO_NO_DIFF<T,TV>& b,const AUTO_NO_DIFF<T,TV>& c)
        :x(a.x,b.x,c.x)
    {
    }

    static AUTO_NO_DIFF From_Var(TV v)
    {return AUTO_NO_DIFF(v);}

    static AUTO_NO_DIFF From_Const(T a)
    {return AUTO_NO_DIFF(a);}

    AUTO_NO_DIFF<T,TV> operator()(int i) const
    {return AUTO_NO_DIFF<T,TV>(x(i));}

    void Set_Entry(int i,const AUTO_NO_DIFF<T,TV>& a)
    {x(i)=a.x;}

    AUTO_NO_DIFF operator-() const
    {return AUTO_NO_DIFF(-x);}

    AUTO_NO_DIFF operator+() const
    {return *this;}

    AUTO_NO_DIFF operator+(const AUTO_NO_DIFF& a) const
    {return AUTO_NO_DIFF(x+a.x);}

    AUTO_NO_DIFF operator-(const AUTO_NO_DIFF& a) const
    {return AUTO_NO_DIFF(x-a.x);}

    AUTO_NO_DIFF operator*(const AUTO_NO_DIFF<T,TV>& a) const
    {return AUTO_NO_DIFF(x*a.x);}

    AUTO_NO_DIFF operator/(const AUTO_NO_DIFF<T,TV>& a) const
    {return *this*((T)1/a);}

    AUTO_NO_DIFF operator+(TV a) const
    {return AUTO_NO_DIFF(x+a);}

    AUTO_NO_DIFF operator-(TV a) const
    {return AUTO_NO_DIFF(x-a);}

    AUTO_NO_DIFF operator*(T a) const
    {return AUTO_NO_DIFF(x*a);}

    AUTO_NO_DIFF operator/(T a) const
    {return AUTO_NO_DIFF(x/a);}

    AUTO_NO_DIFF<T,TV> Dot(TV v) const
    {return AUTO_NO_DIFF<T,TV>(x.Dot(v));}

    AUTO_NO_DIFF<T,TV> Dot(const AUTO_NO_DIFF& a) const
    {return AUTO_NO_DIFF<T,TV>(x.Dot(a.x));}

    AUTO_NO_DIFF<T,TV> Magnitude_Squared() const
    {return Dot(*this);}

    AUTO_NO_DIFF<T,TV> Magnitude() const;
    
    AUTO_NO_DIFF<T,TV> Sum() const
    {return AUTO_NO_DIFF<T,TV>(x.Sum());}
};

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,0> > Magnitude_Helper(const AUTO_NO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_NO_DIFF<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,1> > Magnitude_Helper(const AUTO_NO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return abs(a(0));}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,2> > Magnitude_Helper(const AUTO_NO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return hypot(a(0),a(1));}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,3> > Magnitude_Helper(const AUTO_NO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return hypot(a(0),a(1),a(2));}

template<class TV>
inline AUTO_NO_DIFF<typename TV::SCALAR,TV> AUTO_NO_DIFF<TV,TV>::Magnitude() const {return Magnitude_Helper(*this);}

template<class TV>
inline AUTO_NO_DIFF<TV,TV> operator*(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,TV v)
{return AUTO_NO_DIFF<TV,TV>(v*a.x);}

template<class TV>
inline AUTO_NO_DIFF<TV,TV> operator*(TV v,const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a)
{return a*v;}

template<class TV>
inline AUTO_NO_DIFF<TV,TV> operator*(const AUTO_NO_DIFF<typename TV::SCALAR,TV>& a,const AUTO_NO_DIFF<TV,TV>& v)
{return v*a;}

template<class TV>
inline AUTO_NO_DIFF<TV,TV> operator*(typename TV::SCALAR a,const AUTO_NO_DIFF<TV,TV>& v)
{return v*a;}

template<class T>
inline AUTO_NO_DIFF<VECTOR<T,0>,VECTOR<T,0> > abs(const AUTO_NO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_NO_DIFF<VECTOR<T,0>,VECTOR<T,0> >();}

template<class T>
inline AUTO_NO_DIFF<VECTOR<T,1>,VECTOR<T,1> > abs(const AUTO_NO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return AUTO_NO_DIFF<VECTOR<T,1>,VECTOR<T,1> >(abs(a(0)));}

template<class T>
inline AUTO_NO_DIFF<VECTOR<T,2>,VECTOR<T,2> > abs(const AUTO_NO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return AUTO_NO_DIFF<VECTOR<T,2>,VECTOR<T,2> >(abs(a(0)),abs(a(1)));}

template<class T>
inline AUTO_NO_DIFF<VECTOR<T,3>,VECTOR<T,3> > abs(const AUTO_NO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return AUTO_NO_DIFF<VECTOR<T,3>,VECTOR<T,3> >(abs(a(0)),abs(a(1)),abs(a(2)));}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,0> > max(const AUTO_NO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_NO_DIFF<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,1> > max(const AUTO_NO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,2> > max(const AUTO_NO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return max(a(0),a(1));}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,3> > max(const AUTO_NO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return max(a(0),max(a(1),a(2)));}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,0> > min(const AUTO_NO_DIFF<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_NO_DIFF<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,1> > min(const AUTO_NO_DIFF<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,2> > min(const AUTO_NO_DIFF<VECTOR<T,2>,VECTOR<T,2> >& a)
{return min(a(0),a(1));}

template<class T>
inline AUTO_NO_DIFF<T,VECTOR<T,3> > min(const AUTO_NO_DIFF<VECTOR<T,3>,VECTOR<T,3> >& a)
{return min(a(0),min(a(1),a(2)));}

}
#endif
