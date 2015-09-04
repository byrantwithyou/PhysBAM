//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_DIFFERENCE
//#####################################################################
#ifndef __ARRAY_DIFFERENCE__
#define __ARRAY_DIFFERENCE__

#include <Tools/Arrays/ARRAY_SUM.h>
#include <cassert>

#ifdef DIFFERENCE
#undef DIFFERENCE
#endif

namespace PhysBAM{

template<class T_ARRAY0,class T_ARRAY1> class ARRAY_DIFFERENCE;
template<class T_ARRAY0,class T_ARRAY1> struct IS_ARRAY<ARRAY_DIFFERENCE<T_ARRAY0,T_ARRAY1> > {static const bool value=true;};
template<class T_ARRAY0,class T_ARRAY1> struct IS_ARRAY_VIEW<ARRAY_DIFFERENCE<T_ARRAY0,T_ARRAY1> > {static const bool value=true;};

template<class T_ARRAY0,class T_ARRAY1>
class ARRAY_DIFFERENCE:public ARRAY_EXPRESSION<decltype(*(typename T_ARRAY0::ELEMENT*)0+*(typename T_ARRAY1::ELEMENT*)0),ARRAY_DIFFERENCE<T_ARRAY0,T_ARRAY1>,typename T_ARRAY0::INDEX>
{
    typedef typename T_ARRAY0::ELEMENT T0;typedef typename T_ARRAY1::ELEMENT T2;
    typedef typename conditional<IS_ARRAY_VIEW<T_ARRAY0>::value,const T_ARRAY0,const T_ARRAY0&>::type T_ARRAY1_VIEW; // if it's an array view we can copy it, otherwise store a reference
    typedef typename conditional<IS_ARRAY_VIEW<T_ARRAY1>::value,const T_ARRAY1,const T_ARRAY1&>::type T_ARRAY2_VIEW;
    typedef decltype(*(T0*)0-*(T2*)0) T_DIFFERENCE;
public:
    typedef T_DIFFERENCE ELEMENT;typedef typename T_ARRAY0::INDEX INDEX;

    T_ARRAY1_VIEW array0;
    T_ARRAY2_VIEW array1;

    ARRAY_DIFFERENCE(const T_ARRAY0& array0,const T_ARRAY1& array1)
        :array0(array0),array1(array1)
    {}

    INDEX Size() const
    {INDEX size=array0.Size();assert(size==array1.Size());return size;}

    typename DOMAIN_INDEX_TYPE<INDEX>::TYPE Domain_Indices() const
    {typename DOMAIN_INDEX_TYPE<INDEX>::TYPE domain_indices=array0.Domain_Indices();assert(domain_indices==array1.Domain_Indices());return domain_indices;}

    const T_DIFFERENCE operator()(const INDEX i) const
    {return array0(i)-array1(i);}

//#####################################################################
};

template<class T_ARRAY0,class T_ARRAY1,class ENABLE=void> struct ARRAY_DIFFERENCE_VALID {static const bool value=false;};
template<class T_ARRAY0,class T_ARRAY1> struct ARRAY_DIFFERENCE_VALID<T_ARRAY0,T_ARRAY1,typename enable_if<is_same<typename EQUIVALENT_ARRAY<T_ARRAY0>::TYPE,typename EQUIVALENT_ARRAY<T_ARRAY1>::TYPE>::value>::type>
{static const bool value=IS_ARRAY<T_ARRAY0>::value && IS_ARRAY<T_ARRAY1>::value && (!FIXED_SIZE_VECTOR<T_ARRAY0>::value || !FIXED_SIZE_VECTOR<T_ARRAY1>::value);};

template<class T0,class T2,class T_ARRAY0,class T_ARRAY1> typename enable_if<ARRAY_DIFFERENCE_VALID<T_ARRAY0,T_ARRAY1>::value,ARRAY_DIFFERENCE<T_ARRAY0,T_ARRAY1> >::type
operator-(const ARRAY_BASE<T0,T_ARRAY0,typename T_ARRAY1::INDEX>& array0,const ARRAY_BASE<T2,T_ARRAY1,typename T_ARRAY1::INDEX>& array1)
{return ARRAY_DIFFERENCE<T_ARRAY0,T_ARRAY1>(array0.Derived(),array1.Derived());}

//#####################################################################

}
#endif
