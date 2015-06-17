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
    typedef typename IF<is_reference<T_INDICES>::value,const T_INDICES_NO_REFERENCE&,const T_INDICES_NO_REFERENCE>::TYPE CONST_T_INDICES;
    STATIC_ASSERT(is_same<typename T_ARRAY::INDEX,typename T_INDICES_NO_REFERENCE::ELEMENT>::value);
//    STATIC_ASSERT((!is_const<T_INDICES_NO_REFERENCE>::value));
    typedef typename T_ARRAY::ELEMENT T;typedef typename T_INDICES_NO_REFERENCE::INDEX ID;
    typedef ARRAY_BASE<T,INDIRECT_ARRAY<T_ARRAY,T_INDICES>,ID> BASE;
    struct UNUSABLE{};
    typedef typename IF<IS_ARRAY_VIEW<T_ARRAY>::value,T_ARRAY,T_ARRAY&>::TYPE T_ARRAY_VIEW;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T ELEMENT;typedef ID INDEX;
    typedef typename ARRAY_RESULT_TYPE<T_ARRAY>::TYPE RESULT_TYPE;
    typedef typename T_ARRAY::CONST_RESULT_TYPE CONST_RESULT_TYPE;

    T_ARRAY_VIEW array;
    CONST_T_INDICES indices;

//     template<class T_OTHER_ARRAY>
//     INDIRECT_ARRAY(T_OTHER_ARRAY& array,typename ADD_REFERENCE<CONST_T_INDICES>::TYPE indices,typename DISABLE_IF<IS_ARRAY_VIEW<T_OTHER_ARRAY>::value,UNUSABLE>::TYPE unusable=UNUSABLE())
//         :array(array),indices(indices)
//     {
//         STATIC_ASSERT(is_base_of<T_ARRAY,T_OTHER_ARRAY>::value); // avoid grabbing reference to temporary
//     }

//     template<class T_OTHER_ARRAY>
//     INDIRECT_ARRAY(T_OTHER_ARRAY array,typename ADD_REFERENCE<CONST_T_INDICES>::TYPE indices,typename ENABLE_IF<IS_ARRAY_VIEW<T_OTHER_ARRAY>::value,UNUSABLE>::TYPE unusable=UNUSABLE())
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

    typename IF<is_const<T_ARRAY>::value,const T*,T*>::TYPE Get_Array_Pointer()
    {return array.Get_Array_Pointer()+Offset_If_Contiguous(indices);}

    const T* Get_Array_Pointer() const
    {return array.Get_Array_Pointer()+Offset_If_Contiguous(indices);}

    bool Using_Externally_Allocated_Pointer()
    {return array.Using_Externally_Allocated_Pointer();}

private:
    static ID Offset_If_Contiguous(const IDENTITY_ARRAY<ID>& indices) // for contiguous indices, we can extract an offset
    {return ID();}

    template<class T_INDICES_2>
    static ID Offset_If_Contiguous(const ARRAY_PLUS_SCALAR<int,T_INDICES_2>& indices)
    {return indices.c+Offset_If_Contiguous(indices.array);}

    template<class RW> void Read(std::istream& input)
    {PHYSBAM_NOT_IMPLEMENTED();}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,Size());Write_Binary_Array<RW>(output,Get_Array_Pointer(),Size());}
//#####################################################################
};
}
#endif
