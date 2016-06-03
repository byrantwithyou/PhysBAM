//#####################################################################
// Copyright 2007-2016, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_EXPRESSION
//#####################################################################
#ifndef __ARRAY_EXPRESSION__
#define __ARRAY_EXPRESSION__

#include <Tools/Arrays/ARRAY_BASE.h>
#include <Tools/Arrays/SIMPLE_ITERATOR.h>
#include <Tools/Math_Tools/Inverse.h>
#include <cassert>
namespace PhysBAM{

template<class OP,class ID> class ARRAY_EXPRESSION;
template<class OP,class ID> struct IS_ARRAY<ARRAY_EXPRESSION<OP,ID> > {static const bool value=true;};
template<class OP,class ID> struct IS_ARRAY_VIEW<ARRAY_EXPRESSION<OP,ID> > {static const bool value=true;};

template<class A> using SAFE_ARRAY_HOLDER=typename conditional<IS_ARRAY_VIEW<A>::value,const A,const A&>::type;

template<class OP,class ID>
class ARRAY_EXPRESSION:public ARRAY_BASE<decltype((*(OP*)0)(ID())),ARRAY_EXPRESSION<OP,ID>,ID>
{
public:
    typedef ID INDEX;
    typedef decltype((*(OP*)0)(INDEX())) ELEMENT;
    typedef const ELEMENT CONST_RESULT_TYPE;
    typedef const ELEMENT RESULT_TYPE;
    typedef SIMPLE_ITERATOR<ARRAY_EXPRESSION> iterator;
    typedef SIMPLE_ITERATOR<const ARRAY_EXPRESSION> const_iterator;

    OP op;
    INDEX size;

    ARRAY_EXPRESSION(const OP& op,const INDEX& size)
        :op(op),size(size)
    {}

    INDEX Size() const
    {return size;}

    INDEX Domain_Indices() const
    {return size;}

    const ELEMENT operator()(const INDEX i) const
    {return op(i);}

    SIMPLE_ITERATOR<ARRAY_EXPRESSION> begin()
    {return SIMPLE_ITERATOR<ARRAY_EXPRESSION>(*this,0);}

    SIMPLE_ITERATOR<const ARRAY_EXPRESSION> begin() const
    {return SIMPLE_ITERATOR<const ARRAY_EXPRESSION>(*this,0);}

    SIMPLE_ITERATOR<ARRAY_EXPRESSION> end()
    {return SIMPLE_ITERATOR<ARRAY_EXPRESSION>(*this,size);}

    SIMPLE_ITERATOR<const ARRAY_EXPRESSION> end() const
    {return SIMPLE_ITERATOR<const ARRAY_EXPRESSION>(*this,size);}

//#####################################################################
};

template<class OP,class ID> ARRAY_EXPRESSION<OP,ID>
Make_Array_Expression(const OP& op,const ID& size)
{return ARRAY_EXPRESSION<OP,ID>(op,size);}

template<class T0,class T1,class T_ARRAY0,class T_ARRAY1,class ID,class OP> auto
Array_Expression_Helper(const ARRAY_BASE<T0,T_ARRAY0,ID>& array0,const ARRAY_BASE<T1,T_ARRAY1,ID>& array1,const OP& op,
    typename enable_if<!(FIXED_SIZE_VECTOR<T_ARRAY0>::value || FIXED_SIZE_VECTOR<T_ARRAY1>::value),int>::type)
{
    std::tuple<SAFE_ARRAY_HOLDER<T_ARRAY0>,SAFE_ARRAY_HOLDER<T_ARRAY1> > t(array0.Derived(),array1.Derived());
    return Make_Array_Expression([=](ID i){return op(std::get<0>(t)(i),std::get<1>(t)(i));},array0.Size());
}

template<class T0,class T1,class T_ARRAY0,class T_ARRAY1,class OP> auto
Array_Expression_Helper(const ARRAY_BASE<T0,T_ARRAY0,int>& array0,const ARRAY_BASE<T1,T_ARRAY1,int>& array1,const OP& op,
    typename enable_if<FIXED_SIZE_VECTOR<T_ARRAY0>::value || FIXED_SIZE_VECTOR<T_ARRAY1>::value,int>::type)
{
    static const int size=FIXED_SIZE_VECTOR<T_ARRAY0>::size & FIXED_SIZE_VECTOR<T_ARRAY1>::size;
    VECTOR<T0,size> a(array0);
    VECTOR<T1,size> b(array1);
    VECTOR<decltype(op(a(0),b(0))),size> r;
    for(int i=0;i<r.m;i++) r(i)=op(a(i),b(i));
    return r;
}

template<class T,class T_ARRAY,class ID,class OP> auto
Array_Expression_Helper(const ARRAY_BASE<T,T_ARRAY,ID>& array,const OP& op)
{
    std::tuple<SAFE_ARRAY_HOLDER<T_ARRAY> > t(array.Derived());
    return Make_Array_Expression([=](ID i){return op(std::get<0>(t)(i));},array.Size());
}

template<class T0,class T1,class T_ARRAY0,class T_ARRAY1,class ID> auto
operator+(const ARRAY_BASE<T0,T_ARRAY0,ID>& array0,const ARRAY_BASE<T1,T_ARRAY1,ID>& array1)
{return Array_Expression_Helper(array0,array1,[](const T0& a,const T1& b){return a+b;});}

template<class T0,class T1,class T_ARRAY0,class T_ARRAY1,class ID> auto
operator-(const ARRAY_BASE<T0,T_ARRAY0,ID>& array0,const ARRAY_BASE<T1,T_ARRAY1,ID>& array1)
{return Array_Expression_Helper(array0,array1,[](const T0& a,const T1& b){return a-b;});}

template<class T0,class T1,class T_ARRAY0,class T_ARRAY1,class ID> auto
operator*(const ARRAY_BASE<T0,T_ARRAY0,ID>& array0,const ARRAY_BASE<T1,T_ARRAY1,ID>& array1)
{return Array_Expression_Helper(array0,array1,[](const T0& a,const T1& b){return a*b;});}

template<class T0,class T1,class T_ARRAY0,class T_ARRAY1,class ID> auto
operator/(const ARRAY_BASE<T0,T_ARRAY0,ID>& array0,const ARRAY_BASE<T1,T_ARRAY1,ID>& array1)
{return Array_Expression_Helper(array0,array1,[](const T0& a,const T1& b){return a/b;});}

template<class T,class T_ARRAY,class ID> auto
operator*(const typename T_ARRAY::SCALAR& s,const ARRAY_BASE<T,T_ARRAY,ID>& array)
{return Array_Expression_Helper(array,[=](const T& x){return s*x;});}

template<class T,class T_ARRAY,class ID> auto
operator*(const ARRAY_BASE<T,T_ARRAY,ID>& array,const typename T_ARRAY::SCALAR& s)
{return Array_Expression_Helper(array,[=](const T& x){return x*s;});}

template<class T,class T_ARRAY,class ID> auto
operator/(const ARRAY_BASE<T,T_ARRAY,ID>& array,const typename T_ARRAY::SCALAR& s)
{auto I=Inverse(s);return Array_Expression_Helper(array,[=](const T& x){return x*I;});}

template<class T,class T_ARRAY,class ID> auto
operator-(const ARRAY_BASE<T,T_ARRAY,ID>& array)
{return Array_Expression_Helper(array,[=](const T& x){return -x;});}

template<class T,class T_ARRAY,class ID> auto
operator+(const ARRAY_BASE<T,T_ARRAY,ID>& array,const typename T_ARRAY::SCALAR& s)
{return Array_Expression_Helper(array,[=](const T& x){return x+s;});}

template<class T,class T_ARRAY,class ID> auto
operator+(const ARRAY_BASE<T,T_ARRAY,ID>& array,const typename enable_if<!is_same<typename T_ARRAY::SCALAR,T>::value,T>::type& s)
{return Array_Expression_Helper(array,[=](const T& x){return x+s;});}

template<class T,class T_ARRAY,class ID> auto
operator-(const ARRAY_BASE<T,T_ARRAY,ID>& array,const typename T_ARRAY::SCALAR& s)
{return Array_Expression_Helper(array,[=](const T& x){return x-s;});}

template<class T,class T_ARRAY,class ID> auto
operator-(const ARRAY_BASE<T,T_ARRAY,ID>& array,const typename enable_if<!is_same<typename T_ARRAY::SCALAR,T>::value,T>::type& s)
{return Array_Expression_Helper(array,[=](const T& x){return x-s;});}

template<class T,class T_ARRAY,class ID> auto
operator+(const typename T_ARRAY::SCALAR& s,const ARRAY_BASE<T,T_ARRAY,ID>& array)
{return Array_Expression_Helper(array,[=](const T& x){return s+x;});}

template<class T,class T_ARRAY,class ID> auto
operator+(const typename enable_if<!is_same<typename T_ARRAY::SCALAR,T>::value,T>::type& s,const ARRAY_BASE<T,T_ARRAY,ID>& array)
{return Array_Expression_Helper(array,[=](const T& x){return s+x;});}

template<class T,class T_ARRAY,class ID> auto
operator-(const typename T_ARRAY::SCALAR& s,const ARRAY_BASE<T,T_ARRAY,ID>& array)
{return Array_Expression_Helper(array,[=](const T& x){return s-x;});}

template<class T,class T_ARRAY,class ID> auto
operator-(const typename enable_if<!is_same<typename T_ARRAY::SCALAR,T>::value,T>::type& s,const ARRAY_BASE<T,T_ARRAY,ID>& array)
{return Array_Expression_Helper(array,[=](const T& x){return s-x;});}
}
#endif
