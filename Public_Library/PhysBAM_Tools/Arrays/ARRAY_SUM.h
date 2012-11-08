//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_SUM
//#####################################################################
#ifndef __ARRAY_SUM__
#define __ARRAY_SUM__

#include <PhysBAM_Tools/Arrays/ARRAY_BASE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_EXPRESSION.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
#include <cassert>
namespace PhysBAM{

template<class T_ARRAY1,class T_ARRAY2> class ARRAY_SUM;
template<class T_ARRAY1,class T_ARRAY2> struct IS_ARRAY<ARRAY_SUM<T_ARRAY1,T_ARRAY2> > {static const bool value=true;};
template<class T_ARRAY1,class T_ARRAY2> struct IS_ARRAY_VIEW<ARRAY_SUM<T_ARRAY1,T_ARRAY2> > {static const bool value=true;};

template<class T_ARRAY1,class T_ARRAY2>
class ARRAY_SUM:public ARRAY_EXPRESSION<typename SUM<typename T_ARRAY1::ELEMENT,typename T_ARRAY2::ELEMENT>::TYPE,ARRAY_SUM<T_ARRAY1,T_ARRAY2>,typename T_ARRAY1::INDEX>
{
    typedef typename T_ARRAY1::ELEMENT T1;typedef typename T_ARRAY2::ELEMENT T2;
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY1>::value,const T_ARRAY1,const T_ARRAY1&>::TYPE T_ARRAY1_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY2>::value,const T_ARRAY2,const T_ARRAY2&>::TYPE T_ARRAY2_VIEW;
    typedef typename SUM<T1,T2>::TYPE T_SUM;
public:
    typedef T_SUM ELEMENT;typedef typename T_ARRAY1::INDEX INDEX;

    T_ARRAY1_VIEW array1;
    T_ARRAY2_VIEW array2;

    ARRAY_SUM(const T_ARRAY1& array1,const T_ARRAY2& array2)
        :array1(array1),array2(array2)
    {}

    INDEX Size() const
    {INDEX size=array1.Size();assert(size==array2.Size());return size;}

    typename DOMAIN_INDEX_TYPE<INDEX>::TYPE Domain_Indices() const
    {typename DOMAIN_INDEX_TYPE<INDEX>::TYPE domain_indices=array1.Domain_Indices();assert(domain_indices==array2.Domain_Indices());return domain_indices;}

    const T_SUM operator()(const INDEX i) const
    {return array1(i)+array2(i);}

//#####################################################################
};

template<class T_ARRAY1,class T_ARRAY2,class ENABLE=void> struct ARRAY_SUM_VALID {static const bool value=false;};
template<class T_ARRAY1,class T_ARRAY2> struct ARRAY_SUM_VALID<T_ARRAY1,T_ARRAY2,typename ENABLE_IF<IS_SAME<typename EQUIVALENT_ARRAY<T_ARRAY1>::TYPE,typename EQUIVALENT_ARRAY<T_ARRAY2>::TYPE>::value>::TYPE>
{static const bool value=IS_ARRAY<T_ARRAY1>::value && IS_ARRAY<T_ARRAY2>::value && (!FIXED_SIZE_VECTOR<T_ARRAY1>::value || !FIXED_SIZE_VECTOR<T_ARRAY2>::value);};

template<class T1,class T2,class T_ARRAY1,class T_ARRAY2> typename ENABLE_IF<ARRAY_SUM_VALID<T_ARRAY1,T_ARRAY2>::value,ARRAY_SUM<T_ARRAY1,T_ARRAY2> >::TYPE
operator+(const ARRAY_BASE<T1,T_ARRAY1,typename T_ARRAY2::INDEX>& array1,const ARRAY_BASE<T2,T_ARRAY2,typename T_ARRAY2::INDEX>& array2)
{return ARRAY_SUM<T_ARRAY1,T_ARRAY2>(array1.Derived(),array2.Derived());}

//#####################################################################

template<class T_ARRAY1,class T_ARRAY2> struct SUM<T_ARRAY1,T_ARRAY2,typename ENABLE_IF<ARRAY_SUM_VALID<T_ARRAY1,T_ARRAY2>::value>::TYPE>
{typedef ARRAY_SUM<T_ARRAY1,T_ARRAY2> TYPE;};

//#####################################################################
}
#endif
