//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_NEGATION
//#####################################################################
#ifndef __ARRAY_NEGATION__
#define __ARRAY_NEGATION__

#include <Tools/Arrays/ARRAY_EXPRESSION.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Vectors/ARITHMETIC_POLICY.h>
namespace PhysBAM{

template<class T_ARRAY> class ARRAY_NEGATION;
template<class T_ARRAY> struct IS_ARRAY<ARRAY_NEGATION<T_ARRAY> > {static const bool value=true;};
template<class T_ARRAY> struct IS_ARRAY_VIEW<ARRAY_NEGATION<T_ARRAY> > {static const bool value=true;};

template<class T_ARRAY>
class ARRAY_NEGATION:public ARRAY_EXPRESSION<typename T_ARRAY::ELEMENT,ARRAY_NEGATION<T_ARRAY>,typename T_ARRAY::INDEX>
{
    typedef typename T_ARRAY::ELEMENT T;
    typedef typename conditional<IS_ARRAY_VIEW<T_ARRAY>::value,const T_ARRAY,const T_ARRAY&>::type T_ARRAY_VIEW; // if it's an array view we can copy it, otherwise store a reference
public:
    typedef typename NEGATION<T>::TYPE ELEMENT;typedef typename T_ARRAY::INDEX INDEX;

    T_ARRAY_VIEW array;

    explicit ARRAY_NEGATION(const T_ARRAY& array)
        :array(array)
    {}

    INDEX Size() const
    {return array.Size();}

    typename DOMAIN_INDEX_TYPE<INDEX>::TYPE Domain_Indices() const
    {return array.Domain_Indices();}

    const ELEMENT operator()(const INDEX i) const
    {return -array(i);}

//#####################################################################
};

template<class T_ARRAY1,class ENABLE=void> struct ARRAY_NEGATION_VALID {static const bool value=false;};
template<class T_ARRAY1> struct ARRAY_NEGATION_VALID<T_ARRAY1,typename FIRST<void,typename NEGATION<typename T_ARRAY1::ELEMENT>::TYPE>::TYPE>
{static const bool value=!FIXED_SIZE_VECTOR<T_ARRAY1>::value && IS_ARRAY<T_ARRAY1>::value;};

template<class T,class T_ARRAY> typename enable_if<ARRAY_NEGATION_VALID<T_ARRAY>::value,ARRAY_NEGATION<T_ARRAY> >::type
operator-(const ARRAY_BASE<T,T_ARRAY,typename T_ARRAY::INDEX>& array)
{return ARRAY_NEGATION<T_ARRAY>(array.Derived());}

//#####################################################################

template<class T_ARRAY> struct NEGATION<T_ARRAY,typename enable_if<ARRAY_NEGATION_VALID<T_ARRAY>::value>::type>
{typedef ARRAY_NEGATION<T_ARRAY> TYPE;};

//#####################################################################

}
#endif
