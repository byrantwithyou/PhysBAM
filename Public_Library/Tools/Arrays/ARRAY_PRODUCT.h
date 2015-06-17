//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_PRODUCT
//#####################################################################
#ifndef __ARRAY_PRODUCT__
#define __ARRAY_PRODUCT__

#include <Tools/Arrays/ARRAY_EXPRESSION.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Vectors/ARITHMETIC_POLICY.h>
#include <cassert>
namespace PhysBAM{

template<class T_ARRAY0,class T_ARRAY1> class ARRAY_PRODUCT;
template<class T_ARRAY0,class T_ARRAY1> struct IS_ARRAY<ARRAY_PRODUCT<T_ARRAY0,T_ARRAY1> > {static const bool value=true;};
template<class T_ARRAY0,class T_ARRAY1> struct IS_ARRAY_VIEW<ARRAY_PRODUCT<T_ARRAY0,T_ARRAY1> > {static const bool value=true;};

template<class T_ARRAY0,class T_ARRAY1>
class ARRAY_PRODUCT:public ARRAY_EXPRESSION<typename PRODUCT<typename T_ARRAY0::ELEMENT,typename T_ARRAY1::ELEMENT>::TYPE,ARRAY_PRODUCT<T_ARRAY0,T_ARRAY1>,typename T_ARRAY0::INDEX>
{
    typedef typename T_ARRAY0::ELEMENT T0;typedef typename T_ARRAY1::ELEMENT T2;
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY0>::value,const T_ARRAY0,const T_ARRAY0&>::TYPE T_ARRAY1_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY1>::value,const T_ARRAY1,const T_ARRAY1&>::TYPE T_ARRAY2_VIEW;
    typedef typename PRODUCT<T0,T2>::TYPE T_PRODUCT;
public:
    typedef T_PRODUCT ELEMENT;typedef typename T_ARRAY0::INDEX INDEX;

    T_ARRAY1_VIEW array0;
    T_ARRAY2_VIEW array1;

    ARRAY_PRODUCT(const T_ARRAY0& array0,const T_ARRAY1& array1)
        :array0(array0),array1(array1)
    {}

    INDEX Size() const
    {INDEX size=array0.Size();assert(size==array1.Size());return size;}

    typename DOMAIN_INDEX_TYPE<INDEX>::TYPE Domain_Indices() const
    {typename DOMAIN_INDEX_TYPE<INDEX>::TYPE domain_indices=array0.Domain_Indices();assert(domain_indices==array1.Domain_Indices());return domain_indices;}

    const T_PRODUCT operator()(const INDEX i) const
    {return array0(i)*array1(i);}

//#####################################################################
};

template<class T_ARRAY0,class T_ARRAY1,class ENABLE=void> struct ARRAY_PRODUCT_VALID {static const bool value=false;};
template<class T_ARRAY0,class T_ARRAY1> struct ARRAY_PRODUCT_VALID<T_ARRAY0,T_ARRAY1,typename ENABLE_IF<is_same<typename T_ARRAY0::ELEMENT,typename T_ARRAY1::ELEMENT>::value>::TYPE>
{static const bool value=IS_ARRAY<T_ARRAY0>::value && IS_ARRAY<T_ARRAY1>::value && (!FIXED_SIZE_VECTOR<T_ARRAY0>::value || !FIXED_SIZE_VECTOR<T_ARRAY1>::value);};

template<class T0,class T2,class T_ARRAY0,class T_ARRAY1> typename ENABLE_IF<ARRAY_PRODUCT_VALID<T_ARRAY0,T_ARRAY1>::value,ARRAY_PRODUCT<T_ARRAY0,T_ARRAY1> >::TYPE
operator*(const ARRAY_BASE<T0,T_ARRAY0,typename T_ARRAY0::INDEX>& array0,const ARRAY_BASE<T2,T_ARRAY1,typename T_ARRAY0::INDEX>& array1)
{return ARRAY_PRODUCT<T_ARRAY0,T_ARRAY1>(array0.Derived(),array1.Derived());}

//#####################################################################

template<class T_ARRAY0,class T_ARRAY1> struct PRODUCT<T_ARRAY0,T_ARRAY1,typename ENABLE_IF<ARRAY_PRODUCT_VALID<T_ARRAY0,T_ARRAY1>::value>::TYPE>
{typedef ARRAY_PRODUCT<T_ARRAY0,T_ARRAY1> TYPE;};

//#####################################################################

}
#endif
