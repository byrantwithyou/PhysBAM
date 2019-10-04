//#####################################################################
// Copyright 2006, Geoffrey Irving, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function pow
//#####################################################################
//
// pow is slow
//
//#####################################################################
#ifndef __pow__
#define __pow__

#include <Core/Math_Tools/cbrt.h>
#include <cmath>
#include <type_traits>
namespace PhysBAM{

using ::std::sqrt;
using ::std::pow;
using ::std::enable_if_t;

template<class T,int n,int d=1,class VOID=void> struct POW_HELPER;

template<class T,int n>
struct POW_HELPER<T,n,2>
{
    static T pow(T a){return sqrt(POW_HELPER<T,n>::pow(a));}
};

template<class T,int n>
struct POW_HELPER<T,n,3>
{
    static T pow(T a){return cbrt(POW_HELPER<T,n>::pow(a));}
};

template<class T>
struct POW_HELPER<T,0>
{
    static T pow(T a){return 1;}
};

template<class T>
struct POW_HELPER<T,1>
{
    static T pow(T a){return a;}
};

template<class T,int n>
struct POW_HELPER<T,n,1,enable_if_t<(n<0)> >
{
    static T pow(T a){return 1/POW_HELPER<T,-n>::pow(a);}
};

template<class T,int n>
struct POW_HELPER<T,n,1,enable_if_t<(n>0 && n%2==0)> >
{
    static T pow(T a){T x=POW_HELPER<T,n/2>::pow(a);return x*x;}
};

template<class T,int n>
struct POW_HELPER<T,n,1,enable_if_t<(n>0 && n%2==1)> >
{
    static T pow(T a){T x=POW_HELPER<T,n/2>::pow(a);return x*x*a;}
};

template<int n,class T> inline T pow(T a)
{
    return POW_HELPER<T,n,1>::pow(a);
}

template<int n,int d,class T> inline T pow(T a)
{
    return POW_HELPER<T,n,d>::pow(a);
}

}
#endif
