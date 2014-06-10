//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_HESS
//##################################################################### 
#ifndef __AUTO_HESS__
#define __AUTO_HESS__

#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

template<class T,class TV> struct AUTO_HESS;

template<class TV>
struct AUTO_HESS<typename TV::SCALAR,TV>
{
    typedef typename TV::SCALAR T;

    T x;
    TV dx;
    SYMMETRIC_MATRIX<T,TV::m> ddx;

    AUTO_HESS():
        x(T())
    {}

    AUTO_HESS(T x,TV dx):
        x(x),dx(dx)
    {}

    AUTO_HESS(T x,TV dx,const SYMMETRIC_MATRIX<T,TV::m>& ddx):
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
    {return AUTO_HESS(x*a.x,a.x*dx+x*a.dx,a.x*ddx+x*a.ddx+SYMMETRIC_MATRIX<T,TV::m>::Symmetric_Outer_Product(dx,a.dx));}

    AUTO_HESS operator/(const AUTO_HESS& a) const
    {T z=x/a.x,ax2=sqr(a.x);TV q=dx-z*a.dx;return AUTO_HESS(z,q/a.x,ddx/a.x-z*a.ddx/a.x-SYMMETRIC_MATRIX<T,TV::m>::Symmetric_Outer_Product(q,a.dx/ax2));}

    AUTO_HESS operator+(T a) const
    {return AUTO_HESS(x+a,dx,ddx);}

    AUTO_HESS operator-(T a) const
    {return AUTO_HESS(x-a,dx,ddx);}

    AUTO_HESS operator*(T a) const
    {return AUTO_HESS(x*a,a*dx,a*ddx);}

    AUTO_HESS operator/(T a) const
    {return AUTO_HESS(x/a,dx/a,ddx/a);}
};

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> operator+(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return AUTO_HESS<typename TV::SCALAR,TV>(c+a.x,a.dx,a.ddx);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> operator-(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return AUTO_HESS<typename TV::SCALAR,TV>(c-a.x,-a.dx,-a.ddx);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> operator*(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return AUTO_HESS<typename TV::SCALAR,TV>(c*a.x,c*a.dx,c*a.ddx);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> operator/(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR z=c/a.x,w=z/a.x;return AUTO_HESS<typename TV::SCALAR,TV>(z,-w*a.dx,2*w/a.x*SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx)-w*a.ddx);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> sqrt(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=sqrt(a.x);TV t=a.dx/(2*s);return AUTO_HESS<typename TV::SCALAR,TV>(s,t,a.ddx/(2*s)-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(t)/s);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> sqr(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return AUTO_HESS<typename TV::SCALAR,TV>(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx)*2);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> max(const AUTO_HESS<typename TV::SCALAR,TV>& a,const AUTO_HESS<typename TV::SCALAR,TV>& b)
{return a.x>b.x?a:b;}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> min(const AUTO_HESS<typename TV::SCALAR,TV>& a,const AUTO_HESS<typename TV::SCALAR,TV>& b)
{return a.x>b.x?b:a;}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> log(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{TV z=a.dx/a.x;return AUTO_HESS<typename TV::SCALAR,TV>(::std::log(a.x),z,a.ddx/a.x-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(z));}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> exp(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=exp(a.x);return AUTO_HESS<typename TV::SCALAR,TV>(s,s*a.dx,s*(a.ddx+SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx)));}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> sin(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return AUTO_HESS<typename TV::SCALAR,TV>(s,c*a.dx,c*a.ddx-s*SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx));}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> cos(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR s=sin(a.x),c=cos(a.x);return AUTO_HESS<typename TV::SCALAR,TV>(c,-s*a.dx,-s*a.ddx-c*SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx));}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> tan(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{typename TV::SCALAR t=::std::tan(a.x),s=1+t*t;return AUTO_HESS<typename TV::SCALAR,TV>(t,s*a.dx,s*a.ddx+2*s*t*SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx));}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> abs(const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return a.x>=0?a:-a;}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(const AUTO_HESS<typename TV::SCALAR,TV>& a,const AUTO_HESS<typename TV::SCALAR,TV>& b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return AUTO_HESS<typename TV::SCALAR,TV>(c,d*a.dx+e*b.dx,d*a.ddx+e*b.ddx+SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(d*b.dx-e*a.dx)/c);
}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(const AUTO_HESS<typename TV::SCALAR,TV>& a,typename TV::SCALAR b)
{
    typename TV::SCALAR c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return AUTO_HESS<typename TV::SCALAR,TV>(c,d*a.dx,d*a.ddx+SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(e*a.dx)/c);
}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return hypot(a,b);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(const AUTO_HESS<typename TV::SCALAR,TV>& a,const AUTO_HESS<typename TV::SCALAR,TV>& b,const AUTO_HESS<typename TV::SCALAR,TV>& c)
{
    typedef SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> SM;
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    SM ab=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(aa*b.dx-bb*a.dx);
    SM bc=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(bb*c.dx-cc*b.dx);
    SM ca=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(cc*a.dx-aa*c.dx);
    return AUTO_HESS<typename TV::SCALAR,TV>(s,aa*a.dx+bb*b.dx+cc*c.dx,aa*a.ddx+bb*b.ddx+cc*c.ddx+(ab+bc+ca)/s);
}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(const AUTO_HESS<typename TV::SCALAR,TV>& a,const AUTO_HESS<typename TV::SCALAR,TV>& b,typename TV::SCALAR c)
{
    typedef SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> SM;
    typename TV::SCALAR s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    SM ab=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(aa*b.dx-bb*a.dx);
    SM bc=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(b.dx);
    SM ca=SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx);
    return AUTO_HESS<typename TV::SCALAR,TV>(s,aa*a.dx+bb*b.dx,aa*a.ddx+bb*b.ddx+ab/s+(bc+ca)*(sqr(cc)/s));
}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(const AUTO_HESS<typename TV::SCALAR,TV>& a,typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV>& b)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV>& a,const AUTO_HESS<typename TV::SCALAR,TV>& b)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(const AUTO_HESS<typename TV::SCALAR,TV>& a,typename TV::SCALAR b,typename TV::SCALAR c)
{
    typename TV::SCALAR t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return AUTO_HESS<typename TV::SCALAR,TV>(s,d*a.dx,d*a.ddx+SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(a.dx)*e);
}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,const AUTO_HESS<typename TV::SCALAR,TV>& a,typename TV::SCALAR c)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> hypot(typename TV::SCALAR b,typename TV::SCALAR c,const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return hypot(a,b,c);}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> atan2(const AUTO_HESS<typename TV::SCALAR,TV>& y,const AUTO_HESS<typename TV::SCALAR,TV>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    TV f=d*y.dx-e*x.dx;
    return AUTO_HESS<typename TV::SCALAR,TV>(::std::atan2(y.x,x.x),f,d*y.ddx-e*x.ddx-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(f,d*x.dx+e*y.dx));
}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> atan2(const AUTO_HESS<typename TV::SCALAR,TV>& y,typename TV::SCALAR x)
{
    typename TV::SCALAR c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    TV f=d*y.dx;
    return AUTO_HESS<typename TV::SCALAR,TV>(::std::atan2(y.x,x),f,d*y.ddx-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(f,e*y.dx));
}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> atan2(typename TV::SCALAR y,const AUTO_HESS<typename TV::SCALAR,TV>& x)
{
    typename TV::SCALAR c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    TV f=-e*x.dx;
    return AUTO_HESS<typename TV::SCALAR,TV>(::std::atan2(y,x.x),f,-e*x.ddx-SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m>::Symmetric_Outer_Product(f,d*x.dx));
}

template<class TV>
struct AUTO_HESS<TV,TV>
{
    typedef typename TV::SCALAR T;

    TV x;
    MATRIX<T,TV::m> dx;
    VECTOR<SYMMETRIC_MATRIX<T,TV::m>,TV::m> ddx;

    AUTO_HESS():
        x(T())
    {}

    AUTO_HESS(TV x,const MATRIX<T,TV::m>& dx):
        x(x),dx(dx)
    {}

    AUTO_HESS(TV x,const MATRIX<T,TV::m>& dx,const VECTOR<SYMMETRIC_MATRIX<T,TV::m>,TV::m>& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    AUTO_HESS(const AUTO_HESS<T,TV>& a)
        :x(a.x),dx(a.dx),ddx(a.ddx)
    {}

    AUTO_HESS(const AUTO_HESS<T,TV>& a,const AUTO_HESS<T,TV>& b)
        :x(a.x,b.x),dx(a.dx,b.dx),ddx(a.ddx,b.ddx)
    {
        dx=dx.Transposed();
    }

    AUTO_HESS(const AUTO_HESS<T,TV>& a,const AUTO_HESS<T,TV>& b,const AUTO_HESS<T,TV>& c)
        :x(a.x,b.x,c.x),dx(a.dx,b.dx,c.dx),ddx(a.ddx,b.ddx,c.ddx)
    {
        dx=dx.Transposed();
    }

    static AUTO_HESS From_Var(TV v)
    {return AUTO_HESS(v,MATRIX<T,TV::m>::Identity_Matrix());}

    AUTO_HESS<T,TV> operator()(int i) const
    {return AUTO_HESS<T,TV>(x(i),dx.Row(i),ddx(i));}

    void Set_Entry(int i,const AUTO_HESS<T,TV>& a)
    {x(i)=a.x;dx.Set_Row(i,a.dx);ddx(i)=a.ddx;}

    AUTO_HESS operator-() const
    {return AUTO_HESS(-x,-dx,-ddx);}

    AUTO_HESS operator+() const
    {return *this;}

    AUTO_HESS operator+(const AUTO_HESS& a) const
    {return AUTO_HESS(x+a.x,dx+a.dx,ddx+a.ddx);}

    AUTO_HESS operator-(const AUTO_HESS& a) const
    {return AUTO_HESS(x-a.x,dx-a.dx,ddx-a.ddx);}

    AUTO_HESS operator*(const AUTO_HESS<T,TV>& a) const
    {AUTO_HESS r(x*a.x,a.x*dx+MATRIX<T,TV::m>::Outer_Product(x,a.dx));
    for(int i=0;i<TV::m;i++)
        r.ddx(i)=ddx(i)*a.x+x(i)*a.ddx+SYMMETRIC_MATRIX<T,TV::m>::Symmetric_Outer_Product(dx.Row(i),a.dx);
    return r;}

    AUTO_HESS operator/(const AUTO_HESS<T,TV>& a) const
    {return *this*((T)1/a);}

    AUTO_HESS operator+(TV a) const
    {return AUTO_HESS(x+a,dx,ddx);}

    AUTO_HESS operator-(TV a) const
    {return AUTO_HESS(x-a,dx,ddx);}

    AUTO_HESS operator*(T a) const
    {AUTO_HESS r(x*a,a*dx);for(int i=0;i<TV::m;i++) r.ddx(i)=ddx(i)*a;return r;}

    AUTO_HESS operator/(T a) const
    {AUTO_HESS r(x/a,dx/a);for(int i=0;i<TV::m;i++) r.ddx(i)=ddx(i)/a;return r;}

    AUTO_HESS<T,TV> Dot(TV v) const
    {
        AUTO_HESS<T,TV> r(x.Dot(v),dx*v);
        for(int i=0;i<TV::m;i++) r.ddx+=ddx(i)*v(i);
        return r;
    }

    AUTO_HESS<T,TV> Dot(const AUTO_HESS& a) const
    {
        AUTO_HESS<T,TV> r(x.Dot(a.x),dx*a.x+a.dx*x);
        for(int i=0;i<TV::m;i++)
            r.ddx+=ddx(i)*a.x(i)+a.ddx(i)*x(i)+
                SYMMETRIC_MATRIX<T,TV::m>::Symmetric_Outer_Product(dx.Row(i),a.dx.Row(i));
        return r;
    }

    AUTO_HESS<T,TV> Magnitude_Squared() const
    {return Dot(*this);}

    AUTO_HESS<T,TV> Magnitude() const;
};

template<class T>
inline AUTO_HESS<T,VECTOR<T,0> > Magnitude_Helper(const AUTO_HESS<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_HESS<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_HESS<T,VECTOR<T,1> > Magnitude_Helper(const AUTO_HESS<VECTOR<T,1>,VECTOR<T,1> >& a)
{return abs(a(0));}

template<class T>
inline AUTO_HESS<T,VECTOR<T,2> > Magnitude_Helper(const AUTO_HESS<VECTOR<T,2>,VECTOR<T,2> >& a)
{return hypot(a(0),a(1));}

template<class T>
inline AUTO_HESS<T,VECTOR<T,3> > Magnitude_Helper(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,3> >& a)
{return hypot(a(0),a(1),a(2));}

template<class TV>
inline AUTO_HESS<typename TV::SCALAR,TV> AUTO_HESS<TV,TV>::Magnitude() const {return Magnitude_Helper(*this);}

template<class TV>
inline AUTO_HESS<TV,TV> operator*(const AUTO_HESS<typename TV::SCALAR,TV>& a,TV v)
{
    AUTO_HESS<TV,TV> r(v*a.x,MATRIX<typename TV::SCALAR,TV::m>::Outer_Product(v,a.dx));
    for(int i=0;i<TV::m;i++) r.ddx(i)=v(i)*a.ddx;
    return r;
}

template<class TV>
inline AUTO_HESS<TV,TV> operator*(TV v,const AUTO_HESS<typename TV::SCALAR,TV>& a)
{return a*v;}

template<class TV>
inline AUTO_HESS<TV,TV> operator*(const AUTO_HESS<typename TV::SCALAR,TV>& a,const AUTO_HESS<TV,TV>& v)
{return v*a;}

template<class TV>
inline AUTO_HESS<TV,TV> operator*(typename TV::SCALAR a,const AUTO_HESS<TV,TV>& v)
{return v*a;}

template<class T>
inline AUTO_HESS<VECTOR<T,0>,VECTOR<T,0> > abs(const AUTO_HESS<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_HESS<VECTOR<T,0>,VECTOR<T,0> >();}

template<class T>
inline AUTO_HESS<VECTOR<T,1>,VECTOR<T,1> > abs(const AUTO_HESS<VECTOR<T,1>,VECTOR<T,1> >& a)
{return AUTO_HESS<VECTOR<T,1>,VECTOR<T,1> >(abs(a(0)));}

template<class T>
inline AUTO_HESS<VECTOR<T,2>,VECTOR<T,2> > abs(const AUTO_HESS<VECTOR<T,2>,VECTOR<T,2> >& a)
{return AUTO_HESS<VECTOR<T,2>,VECTOR<T,2> >(abs(a(0)),abs(a(1)));}

template<class T>
inline AUTO_HESS<VECTOR<T,3>,VECTOR<T,3> > abs(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,3> >& a)
{return AUTO_HESS<VECTOR<T,3>,VECTOR<T,3> >(abs(a(0)),abs(a(1)),abs(a(2)));}

template<class T>
inline AUTO_HESS<T,VECTOR<T,0> > max(const AUTO_HESS<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_HESS<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_HESS<T,VECTOR<T,1> > max(const AUTO_HESS<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_HESS<T,VECTOR<T,2> > max(const AUTO_HESS<VECTOR<T,2>,VECTOR<T,2> >& a)
{return max(a(0),a(1));}

template<class T>
inline AUTO_HESS<T,VECTOR<T,3> > max(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,3> >& a)
{return max(a(0),max(a(1),a(2)));}

template<class T>
inline AUTO_HESS<T,VECTOR<T,0> > min(const AUTO_HESS<VECTOR<T,0>,VECTOR<T,0> >& a)
{return AUTO_HESS<T,VECTOR<T,0> >();}

template<class T>
inline AUTO_HESS<T,VECTOR<T,1> > min(const AUTO_HESS<VECTOR<T,1>,VECTOR<T,1> >& a)
{return a(0);}

template<class T>
inline AUTO_HESS<T,VECTOR<T,2> > min(const AUTO_HESS<VECTOR<T,2>,VECTOR<T,2> >& a)
{return min(a(0),a(1));}

template<class T>
inline AUTO_HESS<T,VECTOR<T,3> > min(const AUTO_HESS<VECTOR<T,3>,VECTOR<T,3> >& a)
{return min(a(0),min(a(1),a(2)));}

template<class T> inline VECTOR<T,0> Principal_Curvatures(const AUTO_HESS<T,VECTOR<T,1> >& h)
{return VECTOR<T,0>();}

template<class T> inline VECTOR<T,1> Principal_Curvatures(const AUTO_HESS<T,VECTOR<T,2> >& h)
{VECTOR<T,2> t=h.dx.Perpendicular();return VECTOR<T,1>(t.Dot(h.ddx*t));}

template<class T> VECTOR<T,2> Principal_Curvatures(const AUTO_HESS<T,VECTOR<T,3> >& h);
}
#endif
