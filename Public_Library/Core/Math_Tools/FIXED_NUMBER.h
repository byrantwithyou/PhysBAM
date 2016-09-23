//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FIXED_NUMBER
//#####################################################################
#ifndef __FIXED_NUMBER__
#define __FIXED_NUMBER__

#include <Core/Read_Write/READ_WRITE_FUNCTIONS.h>
namespace PhysBAM{

template<class T,int m>
struct FIXED_NUMBER
{
    typedef T SCALAR;
    typedef int HAS_UNTYPED_READ_WRITE;

    bool operator!() const
    {return !m;}

    template<int n>
    FIXED_NUMBER<T,m+n> operator+(const FIXED_NUMBER<T,n>&) const
    {return FIXED_NUMBER<T,m+n>();}

    template<int n>
    FIXED_NUMBER<T,m-n> operator-(const FIXED_NUMBER<T,n>&) const
    {return FIXED_NUMBER<T,m-n>();}

    template<int n>
    FIXED_NUMBER<T,m*n> operator*(const FIXED_NUMBER<T,n>&) const
    {return FIXED_NUMBER<T,m*n>();}

    template<int n>
    FIXED_NUMBER<T,m/n> operator/(const FIXED_NUMBER<T,n>&) const
    {return FIXED_NUMBER<T,m/n>();}

    FIXED_NUMBER<T,-m> operator-() const
    {return FIXED_NUMBER<T,-m>();}

    template<int n>
    bool operator==(const FIXED_NUMBER<T,n>&) const
    {return m==n;}

    FIXED_NUMBER Inverse() const
    {return FIXED_NUMBER();}

    static FIXED_NUMBER<T,1> One()
    {return FIXED_NUMBER<T,1>();}

    template<class RW> void Read(std::istream& input)
    {}

    template<class RW> void Write(std::ostream& output) const
    {}

//#####################################################################
};

template<class T> struct IS_FIXED_NUMBER {static const bool value=false;};
template<class T,int m> struct IS_FIXED_NUMBER<FIXED_NUMBER<T,m> > {static const bool value=true;};

template<class T,class A,int m> inline auto operator*(const A& x,const FIXED_NUMBER<T,m>&) -> typename enable_if<!IS_FIXED_NUMBER<A>::value,decltype(x*(T)m)>::type  {return x*(T)m;}
template<class T,class A,int m> inline auto operator*(const FIXED_NUMBER<T,m>&,const A& x) -> typename enable_if<!IS_FIXED_NUMBER<A>::value,decltype((T)m*x)>::type {return (T)m*x;}
template<class T,class A,int m> inline auto operator/(const A& x,const FIXED_NUMBER<T,m>&) -> typename enable_if<!IS_FIXED_NUMBER<A>::value,decltype(x/(T)m)>::type {return x/(T)m;}
template<class T,class A,int m> inline auto operator/(const FIXED_NUMBER<T,m>&,const A& x) -> typename enable_if<!IS_FIXED_NUMBER<A>::value,decltype((T)m/x)>::type {return (T)m/x;}
template<class T,class A,int m> inline A& operator*=(A& x,const FIXED_NUMBER<T,m>&) {return x*=(T)m;}
template<class T,class A,int m> inline A& operator/=(A& x,const FIXED_NUMBER<T,m>&) {return x/=(T)m;}

template<class T,class A> inline typename enable_if<!IS_FIXED_NUMBER<A>::value,const A&>::type operator*(const A& x,const FIXED_NUMBER<T,1>&) {return x;}
template<class T,class A> inline typename enable_if<!IS_FIXED_NUMBER<A>::value,const A&>::type operator*(const FIXED_NUMBER<T,1>&,const A& x) {return x;}
template<class T,class A> inline typename enable_if<!IS_FIXED_NUMBER<A>::value,const A&>::type operator/(const A& x,const FIXED_NUMBER<T,1>&) {return x;}
template<class T,class A> inline A& operator*=(A& x,const FIXED_NUMBER<T,1>&) {return x;}
template<class T,class A> inline A& operator/=(A& x,const FIXED_NUMBER<T,1>&) {return x;}

template<class T> inline const T& operator+(const T& x,const FIXED_NUMBER<T,0>) {return x;}
template<class T> inline const T& operator+(const FIXED_NUMBER<T,0>,const T& x) {return x;}
template<class T> inline const T& operator-(const T& x,const FIXED_NUMBER<T,0>) {return x;}
template<class T> inline T operator-(const FIXED_NUMBER<T,0>,const T& x) {return -x;}
template<class T> inline FIXED_NUMBER<T,0> operator*(const FIXED_NUMBER<T,0>,const T&) {return FIXED_NUMBER<T,0>();}
template<class T> inline FIXED_NUMBER<T,0> operator/(const FIXED_NUMBER<T,0>,const T&) {return FIXED_NUMBER<T,0>();}
template<class T> inline FIXED_NUMBER<T,0> operator*(const T&,const FIXED_NUMBER<T,0>) {return FIXED_NUMBER<T,0>();}

template<class T,int m> inline T operator+(const T& x,const FIXED_NUMBER<T,m>&) {return x+m;}
template<class T,int m> inline T operator+(const FIXED_NUMBER<T,m>&,const T& x) {return m+x;}
template<class T,int m> inline T operator-(const T& x,const FIXED_NUMBER<T,m>&) {return x-m;}
template<class T,int m> inline T operator-(const FIXED_NUMBER<T,m>&,const T& x) {return m-x;}

template<class T,int m> inline std::ostream& operator<<(std::ostream& output,const FIXED_NUMBER<T,m>&) {return output<<(T)m;}

template<class T,int m> typename enable_if<is_scalar<T>::value>::type Fill_From(T& a,const FIXED_NUMBER<T,m>& b){a=(T)m;}
template<class T,int m> void Fill_From(FIXED_NUMBER<T,m>& a,const FIXED_NUMBER<T,m>& b){a=b;}

template<class T> typename enable_if<is_scalar<T>::value,T>::type Choose(T& a,T& b){a=b;}
template<class T,int m> T Choose(const T& a,const FIXED_NUMBER<T,m>& b);
template<class T,int m> T Choose(const FIXED_NUMBER<T,m>& a,const T& b);
template<class T,int m,int n> T Choose(const FIXED_NUMBER<T,m>& a,const FIXED_NUMBER<T,n>& b);
template<class T,int m> FIXED_NUMBER<T,m> Choose(const FIXED_NUMBER<T,m>& a,const FIXED_NUMBER<T,m>& b);

}
#endif
