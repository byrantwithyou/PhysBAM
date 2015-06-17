//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_SUM
//#####################################################################
#ifndef __ARRAY_SUM__
#define __ARRAY_SUM__

#include <Tools/Arrays/ARRAY_BASE.h>
#include <Tools/Arrays/ARRAY_EXPRESSION.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Vectors/ARITHMETIC_POLICY.h>
#include <cassert>
namespace PhysBAM{

template<class T_ARRAY0,class T_ARRAY1> class ARRAY_SUM;
template<class T_ARRAY0,class T_ARRAY1> struct IS_ARRAY<ARRAY_SUM<T_ARRAY0,T_ARRAY1> > {static const bool value=true;};
template<class T_ARRAY0,class T_ARRAY1> struct IS_ARRAY_VIEW<ARRAY_SUM<T_ARRAY0,T_ARRAY1> > {static const bool value=true;};

template<class T_ARRAY0,class T_ARRAY1>
class ARRAY_SUM:public ARRAY_EXPRESSION<typename SUM<typename T_ARRAY0::ELEMENT,typename T_ARRAY1::ELEMENT>::TYPE,ARRAY_SUM<T_ARRAY0,T_ARRAY1>,typename T_ARRAY0::INDEX>
{
    typedef typename T_ARRAY0::ELEMENT T0;typedef typename T_ARRAY1::ELEMENT T2;
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY0>::value,const T_ARRAY0,const T_ARRAY0&>::TYPE T_ARRAY1_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY1>::value,const T_ARRAY1,const T_ARRAY1&>::TYPE T_ARRAY2_VIEW;
    typedef typename SUM<T0,T2>::TYPE T_SUM;
public:
    typedef T_SUM ELEMENT;typedef typename T_ARRAY0::INDEX INDEX;

    T_ARRAY1_VIEW array0;
    T_ARRAY2_VIEW array1;

    ARRAY_SUM(const T_ARRAY0& array0,const T_ARRAY1& array1)
        :array0(array0),array1(array1)
    {}

    INDEX Size() const
    {INDEX size=array0.Size();assert(size==array1.Size());return size;}

    typename DOMAIN_INDEX_TYPE<INDEX>::TYPE Domain_Indices() const
    {typename DOMAIN_INDEX_TYPE<INDEX>::TYPE domain_indices=array0.Domain_Indices();assert(domain_indices==array1.Domain_Indices());return domain_indices;}

    const T_SUM operator()(const INDEX i) const
    {return array0(i)+array1(i);}

//#####################################################################
};

template<class T_ARRAY0,class T_ARRAY1,class ENABLE=void> struct ARRAY_SUM_VALID {static const bool value=false;};
template<class T_ARRAY0,class T_ARRAY1> struct ARRAY_SUM_VALID<T_ARRAY0,T_ARRAY1,typename ENABLE_IF<is_same<typename EQUIVALENT_ARRAY<T_ARRAY0>::TYPE,typename EQUIVALENT_ARRAY<T_ARRAY1>::TYPE>::value>::TYPE>
{static const bool value=IS_ARRAY<T_ARRAY0>::value && IS_ARRAY<T_ARRAY1>::value && (!FIXED_SIZE_VECTOR<T_ARRAY0>::value || !FIXED_SIZE_VECTOR<T_ARRAY1>::value);};

template<class T0,class T2,class T_ARRAY0,class T_ARRAY1> typename ENABLE_IF<ARRAY_SUM_VALID<T_ARRAY0,T_ARRAY1>::value,ARRAY_SUM<T_ARRAY0,T_ARRAY1> >::TYPE
operator+(const ARRAY_BASE<T0,T_ARRAY0,typename T_ARRAY1::INDEX>& array0,const ARRAY_BASE<T2,T_ARRAY1,typename T_ARRAY1::INDEX>& array1)
{return ARRAY_SUM<T_ARRAY0,T_ARRAY1>(array0.Derived(),array1.Derived());}

//#####################################################################

template<class T_ARRAY0,class T_ARRAY1> struct SUM<T_ARRAY0,T_ARRAY1,typename ENABLE_IF<ARRAY_SUM_VALID<T_ARRAY0,T_ARRAY1>::value>::TYPE>
{typedef ARRAY_SUM<T_ARRAY0,T_ARRAY1> TYPE;};

//#####################################################################
}
#endif
