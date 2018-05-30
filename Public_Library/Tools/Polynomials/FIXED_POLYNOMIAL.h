//#####################################################################
// Copyright 2018, Frank Madrid, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FIXED_POLYNOMIAL
//##################################################################### 
#ifndef __FIXED_POLYNOMIAL__
#define __FIXED_POLYNOMIAL__
#include <Core/Vectors/VECTOR.h>

namespace PhysBAM{

template<class T,int d>
class FIXED_POLYNOMIAL
{
public:
    STATIC_ASSERT(d>=0);
    enum {degree=d};
    VECTOR<T,d+1> c;

    FIXED_POLYNOMIAL<T,d-1> Derivative() const
    {
        FIXED_POLYNOMIAL<T,d-1> p;
        for(int i=0;i<degree;i++) p.c(i)=c(i+1)*(i+1);
        return p;
    }

    template<int n>
    FIXED_POLYNOMIAL<T,d-n> operator/(const FIXED_POLYNOMIAL<T,n>& p) const
    {
        FIXED_POLYNOMIAL<T,d-n> q;
        FIXED_POLYNOMIAL<T,n-1> r;
        Quotient_Remainder(p,q,r);
        return q;
    }
        
    template<int n>
    FIXED_POLYNOMIAL<T,n-1> operator%(const FIXED_POLYNOMIAL<T,n>& p) const
    {
        FIXED_POLYNOMIAL<T,d-n> q;
        FIXED_POLYNOMIAL<T,n-1> r;
        Quotient_Remainder(p,q,r);
        return r;
    }

    FIXED_POLYNOMIAL<T,d> operator*(T x) const
    {return {c*x};}

    FIXED_POLYNOMIAL<T,d> operator/(T x) const
    {return {c/x};}

    FIXED_POLYNOMIAL<T,d> operator-() const
    {return {-c};}

    template<int n>
    FIXED_POLYNOMIAL<T,d>& operator+=(const FIXED_POLYNOMIAL<T,n>& p){return *this=*this+p;} 

    template<int n>
    FIXED_POLYNOMIAL<T,d>& operator-=(const FIXED_POLYNOMIAL<T,n>& p){return *this=*this-p;} 

    bool operator==(const FIXED_POLYNOMIAL<T,d>& p) const
    {return c==p.c;}

    bool operator!=(const FIXED_POLYNOMIAL<T,d>&p) const
    {return c!=p.c;}

    T Value(T x) const
    {
        T y=c(degree);
        for(int i=degree-1;i>=0;i--)
            y=c(i)+x*y;
        return y;
    }

    T operator() (T x) const
    {return Value(x);}

    template<int n>
    void Quotient_Remainder(const FIXED_POLYNOMIAL<T,n>& b,
        FIXED_POLYNOMIAL<T,d-n>& q,FIXED_POLYNOMIAL<T,n-1>& r) const
    {
        assert(b.c(n));
        VECTOR<T,d+1> v=c;
        for(int i=d-n;i>=0;i--){
            T f=v(i+n)/b.c(n);
            q.c(i)=f;
            for(int j=0;j<=n;j++)
                v(i+j)-=f*b.c(j);}
        for(int j=0;j<n;j++)
            r.c(j)=v(j);
    }
};

template<class T,int m,int n>
FIXED_POLYNOMIAL<T,std::max(m,n)> operator+(const FIXED_POLYNOMIAL<T,m>& a,const FIXED_POLYNOMIAL<T,n>& b)
{
    FIXED_POLYNOMIAL<T,std::max(m,n)> p;
    for(int i=0;i<m;i++) p.c(i)=a.c(i);
    for(int i=0;i<n;i++) p.c(i)+=b.c(i);
    return p;
}

template<class T,int m,int n>
FIXED_POLYNOMIAL<T,std::max(m,n)> operator-(const FIXED_POLYNOMIAL<T,m>& a,const FIXED_POLYNOMIAL<T,n>& b)
{
    FIXED_POLYNOMIAL<T,std::max(m,n)> p;
    for(int i=0;i<m;i++) p.c(i)=a.c(i);
    for(int i=0;i<n;i++) p.c(i)-=b.c(i);
    return p;
}

template<class T,int m,int n>
FIXED_POLYNOMIAL<T,m+n> operator*(const FIXED_POLYNOMIAL<T,m>& a,const FIXED_POLYNOMIAL<T,n>& b)
{
    FIXED_POLYNOMIAL<T,m+n> p;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            p.c(i+j)+=a.c(i)*b.c(i);
    return p;
}

template<class T,int d>
FIXED_POLYNOMIAL<T,d> operator*(T x,FIXED_POLYNOMIAL<T,d> p)
{
    return p*x;
}

template<class T,int d>
std::ostream& operator<<(std::ostream& o,const FIXED_POLYNOMIAL<T,d>& p)
{
    return o<<p.c;
}
}
#endif
