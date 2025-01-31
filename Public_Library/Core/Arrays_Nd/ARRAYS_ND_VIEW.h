//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_VIEW
//#####################################################################
#ifndef __ARRAYS_ND_VIEW__
#define __ARRAYS_ND_VIEW__

#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Arrays_Nd/ARRAYS_ND_BASE.h>
#include <Core/Math_Tools/exchange.h>
#include <Core/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Core/Utilities/EXCEPTIONS.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int d>
class ARRAY_VIEW<T,VECTOR<int,d> >:public ARRAYS_ND_BASE<typename remove_const<T>::type,VECTOR<int,d> >
{
    typedef VECTOR<int,d> TV_INT;
    struct UNUSABLE{};
    template<class S> struct COPY_CONST{typedef typename conditional<is_const<T>::value,typename add_const<S>::type,S>::type TYPE;};
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef typename remove_const<T>::type ELEMENT;typedef TV_INT INDEX;
    typedef ARRAYS_ND_BASE<ELEMENT,VECTOR<int,d> > BASE;
    typedef T& RESULT_TYPE;

    using BASE::domain;using BASE::array;
private:
    friend class ARRAY_VIEW<typename conditional<is_const<T>::value,ELEMENT,const ELEMENT>::type,TV_INT>;

    using BASE::Calculate_Acceleration_Constants;

    void Initialize(ELEMENT* raw_data)
    {ARRAY_VIEW<ELEMENT> new_array(raw_data,domain.Size());new_array.Exchange(array);}

public:

    ARRAY_VIEW(const RANGE<TV_INT>& domain=RANGE<TV_INT>::Empty_Box(),ELEMENT* raw_data=0)
        :BASE(domain)
    {Initialize(raw_data);}

    ARRAY_VIEW(const ARRAY_VIEW<ELEMENT,TV_INT>& array_input)
        :BASE(array_input.domain)
    {array.Set((ELEMENT*)array_input.array.Get_Array_Pointer(),domain.Size());} // TODO: why do I need this cast?

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY& array_input,enable_if_t<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && !IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE> unusable=UNUSABLE())
        :BASE(array_input.Domain_Indices())
    {Initialize(array_input.array.Get_Array_Pointer());}

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY array_input,enable_if_t<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE> unusable=UNUSABLE())
        :BASE(array_input.Domain_Indices())
    {Initialize(array_input.array.Get_Array_Pointer());}

    template<class T_ARRAY1>
    ARRAY_VIEW& operator=(const T_ARRAY1& source)
    {assert(Equal_Dimensions(*this,source));array=source.array;return *this;}

    void Exchange(ARRAY_VIEW& other)
    {STATIC_ASSERT(!is_const<T>::value);this->BASE::Exchange(other);} // make ARRAY_VIEW<const T> equivalent to const ARRAY_VIEW<const T>

    static void Exchange(ARRAY_VIEW& array0,ARRAY_VIEW& array1)
    {STATIC_ASSERT(!is_const<T>::value);array0.Exchange(array1);}

    void Set(const RANGE<TV_INT>& domain=RANGE<TV_INT>::Empty_Box(),ELEMENT* raw_data=0)
    {Calculate_Acceleration_Constants(domain);Initialize(raw_data);}

    void Set(const ARRAY_VIEW<ELEMENT,TV_INT>& array)
    {Calculate_Acceleration_Constants(array.domain);Initialize(array.base_pointer);}

    template<class T_ARRAY>
    void Set(T_ARRAY& array,enable_if_t<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && !IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE> unusable=UNUSABLE())
    {Calculate_Acceleration_Constants(array.Domain_Indices());Initialize(array.Get_Array_Pointer());}

    template<class T_ARRAY>
    void Set(T_ARRAY array,enable_if_t<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE> unusable=UNUSABLE())
    {Calculate_Acceleration_Constants(array.Domain_Indices());Initialize(array.Get_Array_Pointer());}

    static bool Same_Array(const ARRAY_VIEW& array0,const ARRAY_VIEW& array1)
    {return array0.Get_Array_Pointer()==array1.Get_Array_Pointer();}

    template<class RW>
    void Read(std::istream& input)
    {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
};
}
#endif
