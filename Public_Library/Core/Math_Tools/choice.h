//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function choice
//#####################################################################
#ifndef __choice__
#define __choice__

#include <cassert>
#include <type_traits>
namespace PhysBAM{
using std::conditional;

template<int i,class T1,class T2,class T3> struct choice_helper;

template<class T> struct choice_helper<0,T,void,void>
{
    T& operator()(const int i,T& a,T& b){assert((unsigned)i<2);T* choices[]={&a,&b};return *choices[i];}
    T& operator()(const int i,T& a,T& b,T& c){assert((unsigned)i<3);T* choices[]={&a,&b,&c};return *choices[i];}
};

template<class T1,class T2> struct choice_helper<0,T1,T2,void> {static T1& helper(T1& a,T2& b){return a;}};
template<class T1,class T2> struct choice_helper<1,T1,T2,void> {static T2& helper(T1& a,T2& b){return b;}};

template<class T1,class T2,class T3> struct choice_helper<0,T1,T2,T3> {static T1& helper(T1& a,T2& b,T3& c){return a;}};
template<class T1,class T2,class T3> struct choice_helper<1,T1,T2,T3> {static T2& helper(T1& a,T2& b,T3& c){return b;}};
template<class T1,class T2,class T3> struct choice_helper<2,T1,T2,T3> {static T3& helper(T1& a,T2& b,T3& c){return c;}};

template<int i,class T1,class T2>
inline typename conditional<i==0,T1&,T2&>::type choice(T1& a,T2& b)
{return choice_helper<i,T1,T2,void>::helper(a,b);}

template<int i,class T1,class T2>
inline typename conditional<i==0,const T1&,const T2&>::type choice(const T1& a,const T2& b)
{return choice_helper<i,const T1,const T2,void>::helper(a,b);}

template<int i,class T1,class T2>
inline typename conditional<i==0,T1&,const T2&>::type choice(T1& a,const T2& b)
{return choice_helper<i,T1,const T2,void>::helper(a,b);}

template<int i,class T1,class T2>
inline typename conditional<i==0,const T1&,T2&>::type choice(const T1& a,T2& b)
{return choice_helper<i,const T1,T2,void>::helper(a,b);}

template<int i,class T1,class T2,class T3>
inline typename conditional<i==0,T1&,typename conditional<i==1,T2&,T3&>::type>::type choice(T1& a,T2& b,T3& c)
{return choice_helper<i,T1,T2,T3>::helper(a,b,c);}
}
#endif
