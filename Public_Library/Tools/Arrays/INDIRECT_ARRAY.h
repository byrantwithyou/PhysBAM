//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Tamar Shinar, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INDIRECT_ARRAY
//#####################################################################
#ifndef __INDIRECT_ARRAY__
#define __INDIRECT_ARRAY__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAY_BASE.h>
#include <Tools/Arrays/SIMPLE_ITERATOR.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T_ARRAY,class T_INDICES> struct IS_ARRAY<INDIRECT_ARRAY<T_ARRAY,T_INDICES> > {static const bool value=true;};
template<class T_ARRAY,class T_INDICES> struct IS_ARRAY_VIEW<INDIRECT_ARRAY<T_ARRAY,T_INDICES> > {static const bool value=true;};

template<class T_INDICES> struct INDIRECT_ARRAY_BASE{};
template<class ID,int d> struct INDIRECT_ARRAY_BASE<VECTOR<ID,d>&>{enum {m=VECTOR<ID,d>::m};};
template<class T_ARRAY,class T_INDICES> struct INDIRECT_ARRAY_BASE<INDIRECT_ARRAY<T_ARRAY,T_INDICES>&>:public INDIRECT_ARRAY_BASE<T_INDICES>{};

template<class T_ARRAY,class T_INDICES> // T_INDICES=ARRAY<int>&
class INDIRECT_ARRAY:public INDIRECT_ARRAY_BASE<T_INDICES>,
                     public ARRAY_BASE<typename T_ARRAY::ELEMENT,INDIRECT_ARRAY<T_ARRAY,T_INDICES>,typename remove_reference<T_INDICES>::type::INDEX>
{
    typedef typename remove_reference<T_INDICES>::type T_INDICES_NO_REFERENCE;
    typedef typename conditional<is_reference<T_INDICES>::value,const T_INDICES_NO_REFERENCE&,const T_INDICES_NO_REFERENCE>::type CONST_T_INDICES;
    STATIC_ASSERT(is_same<typename T_ARRAY::INDEX,typename T_INDICES_NO_REFERENCE::ELEMENT>::value);
//    STATIC_ASSERT((!is_const<T_INDICES_NO_REFERENCE>::value));
    typedef typename T_ARRAY::ELEMENT T;typedef typename T_INDICES_NO_REFERENCE::INDEX ID;
    typedef ARRAY_BASE<T,INDIRECT_ARRAY<T_ARRAY,T_INDICES>,ID> BASE;
    struct UNUSABLE{};
    typedef typename conditional<IS_ARRAY_VIEW<T_ARRAY>::value,T_ARRAY,T_ARRAY&>::type T_ARRAY_VIEW;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T ELEMENT;typedef ID INDEX;
    typedef typename ARRAY_RESULT_TYPE<T_ARRAY>::TYPE RESULT_TYPE;
    typedef typename T_ARRAY::CONST_RESULT_TYPE CONST_RESULT_TYPE;
    typedef SIMPLE_ITERATOR<INDIRECT_ARRAY> iterator;
    typedef SIMPLE_ITERATOR<const INDIRECT_ARRAY> const_iterator;

    T_ARRAY_VIEW array;
    CONST_T_INDICES indices;

//     template<class T_OTHER_ARRAY>
//     INDIRECT_ARRAY(T_OTHER_ARRAY& array,typename ADD_REFERENCE<CONST_T_INDICES>::TYPE indices,typename enable_if<:IS_ARRAY_VIEW<T_OTHER_ARRAY>::value,UNUSABLE>::type unusable=UNUSABLE())
//         :array(array),indices(indices)
//     {
//         STATIC_ASSERT(is_base_of<T_ARRAY,T_OTHER_ARRAY>::value); // avoid grabbing reference to temporary
//     }

//     template<class T_OTHER_ARRAY>
//     INDIRECT_ARRAY(T_OTHER_ARRAY array,typename ADD_REFERENCE<CONST_T_INDICES>::TYPE indices,typename enable_if<IS_ARRAY_VIEW<T_OTHER_ARRAY>::value,UNUSABLE>::type unusable=UNUSABLE())
//         :array(array),indices(indices)
//     {
//     }

    INDIRECT_ARRAY(T_ARRAY_VIEW array,typename add_lvalue_reference<CONST_T_INDICES>::type indices)
        :array(array),indices(indices)
    {
    }

    INDIRECT_ARRAY(const INDIRECT_ARRAY<typename remove_const<T_ARRAY>::type,T_INDICES>& indirect)
        :array(indirect.array),indices(indirect.indices)
    {}

    ID Size() const
    {return indices.Size();}

    RESULT_TYPE operator()(const ID i)
    {return array(indices(i));}

    CONST_RESULT_TYPE operator()(const ID i) const
    {return array(indices(i));}

    INDIRECT_ARRAY& operator=(const INDIRECT_ARRAY& source)
    {return BASE::operator=(source);}

    template<class T_OTHER_ARRAY>
    INDIRECT_ARRAY& operator=(const T_OTHER_ARRAY& source)
    {return BASE::operator=(source);}

    typename conditional<is_const<T_ARRAY>::value,const T*,T*>::type Get_Array_Pointer()
    {return array.Get_Array_Pointer();}

    const T* Get_Array_Pointer() const
    {return array.Get_Array_Pointer();}

    bool Using_Externally_Allocated_Pointer()
    {return array.Using_Externally_Allocated_Pointer();}

    SIMPLE_ITERATOR<INDIRECT_ARRAY> begin()
    {return SIMPLE_ITERATOR<INDIRECT_ARRAY>(*this,0);}

    SIMPLE_ITERATOR<const INDIRECT_ARRAY> begin() const
    {return SIMPLE_ITERATOR<const INDIRECT_ARRAY>(*this,0);}

    SIMPLE_ITERATOR<INDIRECT_ARRAY> end()
    {return SIMPLE_ITERATOR<INDIRECT_ARRAY>(*this,Size());}

    SIMPLE_ITERATOR<const INDIRECT_ARRAY> end() const
    {return SIMPLE_ITERATOR<const INDIRECT_ARRAY>(*this,Size());}

private:
    template<class RW> void Read(std::istream& input)
    {PHYSBAM_NOT_IMPLEMENTED();}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,Size());Write_Binary_Array<RW>(output,Get_Array_Pointer(),Size());}
//#####################################################################
};
}
#endif
