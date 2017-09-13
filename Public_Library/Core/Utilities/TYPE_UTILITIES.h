//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TYPE_UTILITIES
//#####################################################################
#ifndef __TYPE_UTILITIES__
#define __TYPE_UTILITIES__

#include <cstdio>
#include <cstdlib>
#include <type_traits>

namespace PhysBAM{
using std::add_const;
using std::add_lvalue_reference;
using std::add_pointer;
using std::conditional;
using std::enable_if;
using std::is_base_of;
using std::is_class;
using std::is_const;
using std::is_convertible;
using std::is_empty;
using std::is_enum;
using std::is_floating_point;
using std::is_fundamental;
using std::is_integral;
using std::is_pod;
using std::is_pointer;
using std::is_reference;
using std::is_same;
using std::is_scalar;
using std::remove_const;
using std::remove_pointer;
using std::remove_reference;

template<class T1,class T2> struct ASSERT_SAME_HELPER{static const bool value=false;};
template<class T> struct ASSERT_SAME_HELPER<T,T>{static const bool value=true;};

#define STATIC_ASSERT_SAME(T1,T2) static_assert(::PhysBAM::ASSERT_SAME_HELPER<T1,T2>::value,"ASSERT_SAME_HELPER<"#T1","#T2">")

template<class T1,class T2=void,class T3=void,class T4=void> struct FIRST{typedef T1 TYPE;};

template<class T_ARRAY,class ENABLER=void> struct IS_ARRAY {static const bool value=false;};
template<class T_ARRAY> struct IS_ARRAY<const T_ARRAY>:public IS_ARRAY<T_ARRAY>{};

template<class T_ARRAY,class ENABLER=void> struct IS_ARRAY_VIEW {static const bool value=false;};
template<class T_ARRAY> struct IS_ARRAY_VIEW<const T_ARRAY>:public IS_ARRAY_VIEW<T_ARRAY>{};

template<class T,int d> class VECTOR;
template<class TV> struct FIXED_SIZE_VECTOR {static const int value=false;static const int size=-1;};
template<class T,int d> struct FIXED_SIZE_VECTOR<VECTOR<T,d> > {static const int value=true;static const int size=d;};

template<class T> struct HAS_CHEAP_COPY {static const bool value=is_fundamental<T>::value || is_enum<T>::value || IS_ARRAY_VIEW<T>::value;};

template<class T> struct HAS_TRIVIAL_DESTRUCTOR {static const bool value=is_pod<T>::value;};

template<class T,class RW,class ENABLER=void> struct IS_BINARY_IO_SAFE;

template<class T,class SCALAR,class ENABLER=void> struct REPLACE_FLOATING_POINT{};
template<class T,class SCALAR> struct REPLACE_FLOATING_POINT<T,SCALAR,typename enable_if<(is_same<T,float>::value || is_same<T,double>::value) && (is_same<SCALAR,float>::value || is_same<SCALAR,double>::value)>::type>{typedef SCALAR TYPE;};
template<class T,class SCALAR> struct REPLACE_FLOATING_POINT<T,SCALAR,typename enable_if<!(is_same<T,float>::value || is_same<T,double>::value) && is_fundamental<T>::value && (is_same<SCALAR,float>::value || is_same<SCALAR,double>::value)>::type>{typedef T TYPE;};
template<class T,class SCALAR> struct REPLACE_FLOATING_POINT<T,SCALAR,typename enable_if<is_pointer<T>::value>::type> {typedef typename REPLACE_FLOATING_POINT<typename remove_pointer<T>::type,SCALAR>::type* TYPE;};
}
#endif
