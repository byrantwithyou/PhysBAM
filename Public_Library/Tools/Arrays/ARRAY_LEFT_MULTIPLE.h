//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_LEFT_MULTIPLE
//#####################################################################
#ifndef __ARRAY_LEFT_MULTIPLE__
#define __ARRAY_LEFT_MULTIPLE__

#include <Tools/Arrays/ARRAY_EXPRESSION.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Vectors/ARITHMETIC_POLICY.h>
namespace PhysBAM{

template<class T0,class T_ARRAY1> class ARRAY_LEFT_MULTIPLE;
template<class T0,class T_ARRAY1> struct IS_ARRAY<ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> > {static const bool value=true;};
template<class T0,class T_ARRAY1> struct IS_ARRAY_VIEW<ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> > {static const bool value=true;};

template<class T0,class T_ARRAY1>
class ARRAY_LEFT_MULTIPLE:public ARRAY_EXPRESSION<typename PRODUCT<T0,typename T_ARRAY1::ELEMENT>::TYPE,ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1>,typename T_ARRAY1::INDEX>
{
    typedef typename T_ARRAY1::ELEMENT T2;
    typedef typename conditional<HAS_CHEAP_COPY<T0>::value,const T0,const T0&>::type T1_VIEW; // copy if cheap, otherwise store a reference
    typedef typename conditional<IS_ARRAY_VIEW<T_ARRAY1>::value,const T_ARRAY1,const T_ARRAY1&>::type T_ARRAY2_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename PRODUCT<T0,T2>::TYPE T_PRODUCT;
public:
    typedef T_PRODUCT ELEMENT;typedef typename T_ARRAY1::INDEX INDEX;

    T1_VIEW c;
    T_ARRAY2_VIEW array;

    ARRAY_LEFT_MULTIPLE(const T0& c,const T_ARRAY1& array)
        :c(c),array(array)
    {}

    INDEX Size() const
    {return array.Size();}

    typename DOMAIN_INDEX_TYPE<INDEX>::TYPE Domain_Indices() const
    {return array.Domain_Indices();}

    const T_PRODUCT operator()(const INDEX i) const
    {return c*array(i);}

//#####################################################################
};

template<class T0,class T_ARRAY1,class ENABLE=void> struct ARRAY_LEFT_MULTIPLE_VALID {static const bool value=false;};
template<class T0,class T_ARRAY1> struct ARRAY_LEFT_MULTIPLE_VALID<T0,T_ARRAY1,typename FIRST<void,typename PRODUCT<T0,typename T_ARRAY1::ELEMENT>::TYPE>::TYPE>
{static const bool value=!FIXED_SIZE_VECTOR<T_ARRAY1>::value && IS_ARRAY<T_ARRAY1>::value && (is_same<T0,typename T_ARRAY1::ELEMENT>::value || is_scalar<T0>::value);};

template<class T0,class T,class T_ARRAY1,class ID> typename enable_if<ARRAY_LEFT_MULTIPLE_VALID<T0,T_ARRAY1>::value,ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> >::type
operator*(const T0& c,const ARRAY_BASE<T,T_ARRAY1,ID>& array)
{return ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1>(c,array.Derived());}

template<class T0,class T,class T_ARRAY1,class ID> typename enable_if<ARRAY_LEFT_MULTIPLE_VALID<T0,T_ARRAY1>::value,ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> >::type
operator/(const ARRAY_BASE<T,T_ARRAY1,ID>& array,const T0& c)
{STATIC_ASSERT(is_same<T0,float>::value || is_same<T0,double>::value);return ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1>(1/c,array.Derived());}

template<class T0,class T,class T_ARRAY1,class ID> typename enable_if<ARRAY_LEFT_MULTIPLE_VALID<T0,T_ARRAY1>::value,ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> >::type
operator*(const ARRAY_BASE<T,T_ARRAY1,ID>& array,const T0& c)
{STATIC_ASSERT(is_scalar<T0>::value);return ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1>(c,array.Derived());}

//#####################################################################

template<class T0,class T_ARRAY1> struct PRODUCT<T0,T_ARRAY1,typename enable_if<ARRAY_LEFT_MULTIPLE_VALID<T0,T_ARRAY1>::value>::type>
{typedef ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> TYPE;};

template<class T0,class T_ARRAY1> struct PRODUCT<T_ARRAY1,T0,typename enable_if<ARRAY_LEFT_MULTIPLE_VALID<T0,T_ARRAY1>::value>::type>
{typedef ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> TYPE;};

template<class T0,class T_ARRAY1> struct QUOTIENT<T_ARRAY1,T0,typename enable_if<ARRAY_LEFT_MULTIPLE_VALID<T0,T_ARRAY1>::value>::type>
{typedef ARRAY_LEFT_MULTIPLE<T0,T_ARRAY1> TYPE;};

//#####################################################################
}
#endif
