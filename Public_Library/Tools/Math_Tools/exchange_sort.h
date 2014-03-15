//#####################################################################
// Copyright 2002, Robert Bridson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function exchange_sort
//#####################################################################
//
// exchanges the values passed to be in sorted (ascending) order 
//               
//#####################################################################
#ifndef __exchange_sort__
#define __exchange_sort__

#include <Tools/Math_Tools/exchange.h>
namespace PhysBAM{

template<class T>
inline void exchange_sort(T& a,T& b)
{if(a>b) exchange(a,b);}

template<class T>
inline void exchange_sort(T& a,T& b,T& c)
{if(a>b) exchange(a,b);if(b>c) exchange(b,c);if(a>b) exchange(a,b);}

template<class T>
inline void exchange_sort(T& a,T& b,T& c,T& d)
{if(a>b) exchange(a,b);if(c>d) exchange(c,d);if(a>c) exchange(a,c);if(b>d) exchange(b,d);if(b>c) exchange(b,c);}

template<class T,class T_COMPARE>
inline void exchange_sort(T& a,T& b,const T_COMPARE& comp)
{if(comp(b,a)) exchange(a,b);}

template<class T,class T_COMPARE>
inline void exchange_sort(T& a,T& b,T& c,const T_COMPARE& comp)
{if(comp(b,a)) exchange(a,b);if(comp(c,b)) exchange(b,c);if(comp(b,a)) exchange(a,b);}

template<class T,class T_COMPARE>
inline void exchange_sort(T& a,T& b,T& c,T& d,const T_COMPARE& comp)
{if(comp(b,a)) exchange(a,b);if(comp(d,c)) exchange(c,d);if(comp(c,a)) exchange(a,c);if(comp(d,b)) exchange(b,d);if(comp(c,b)) exchange(b,c);}

}
#endif

