//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_VIEW
//#####################################################################
#ifndef __ARRAY_VIEW__
#define __ARRAY_VIEW__

#include <Core/Arrays/ARRAY_BASE.h>
#include <Core/Utilities/EXCEPTIONS.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class ID> struct IS_ARRAY<ARRAY_VIEW<T,ID> > {static const bool value=true;};
template<class T,class ID> struct IS_ARRAY_VIEW<ARRAY_VIEW<T,ID> > {static const bool value=true;};
template<class T,class ID> struct HAS_POINTER<ARRAY_VIEW<T,ID> > {static const bool value=true;};

template<class T,class ID>
class ARRAY_VIEW:public ARRAY_BASE<typename remove_const<T>::type,ARRAY_VIEW<T,ID>,ID>
{
    struct UNUSABLE{};
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef typename remove_const<T>::type ELEMENT;typedef ID INDEX;
    typedef ARRAY_BASE<ELEMENT,ARRAY_VIEW<T,ID>,ID> BASE;
    typedef typename conditional<is_const<T>::value,ELEMENT,UNUSABLE>::type U;
    typedef T& RESULT_TYPE;
    typedef T* iterator; // for stl
    typedef const T* const_iterator; // for stl

    // m and base_pointer inherit constness of T
    T* base_pointer;
    ID m;

    using BASE::Same_Array;

    ARRAY_VIEW()
        :base_pointer(0),m(0)
    {}

    ARRAY_VIEW(T* raw_data,ID m)
        :base_pointer(raw_data),m(m)
    {}

    ARRAY_VIEW(const ARRAY_VIEW<T,ID>&) = default;
    ARRAY_VIEW(ARRAY_VIEW<T,ID>&& a)
        :base_pointer(a.base_pointer),m(a.m)
    {a.base_pointer=0;}

    ARRAY_VIEW(const ARRAY_VIEW<U,ID>& array)
        :base_pointer(array.base_pointer),m(array.m)
    {}

    template<class T_ARRAY>
    ARRAY_VIEW(ARRAY_BASE<ELEMENT,T_ARRAY,ID>& array)
        :base_pointer(array.Derived().Get_Array_Pointer()),m(array.Size())
    {}

    template<class T_ARRAY>
    ARRAY_VIEW(const ARRAY_BASE<ELEMENT,T_ARRAY,ID>& array)
        :base_pointer(array.Derived().Get_Array_Pointer()),m(array.Size())
    {}

    ~ARRAY_VIEW() = default;

    ID Size() const
    {return m;}

    T& operator()(ID i)
    {assert((unsigned)Value(i)<(unsigned)Value(m));return base_pointer[Value(i)];}

    const T& operator()(ID i) const
    {assert((unsigned)Value(i)<(unsigned)Value(m));return base_pointer[Value(i)];}

    bool Valid_Index(ID i) const
    {return (unsigned)Value(i)<(unsigned)Value(m);}

    ARRAY_VIEW& operator=(const ARRAY_VIEW& source)
    {return BASE::operator=(source);}

    ARRAY_VIEW& operator=(ARRAY_VIEW&& source) = delete;

    template<class T_ARRAY1>
    ARRAY_VIEW& operator=(const ARRAY_BASE<ELEMENT,T_ARRAY1,ID>& source)
    {return BASE::operator=(source);}

    T* Get_Array_Pointer()
    {return base_pointer;}

    const T* Get_Array_Pointer() const
    {return base_pointer;}

    T* begin() // for stl
    {return Get_Array_Pointer();}

    const T* begin() const // for stl
    {return Get_Array_Pointer();}

    T* end() // for stl
    {return Get_Array_Pointer()+m;}

    const T* end() const // for stl
    {return Get_Array_Pointer()+m;}

    static bool Same_Array(const ARRAY_VIEW& array0,const ARRAY_VIEW& array1)
    {return array0.Get_Array_Pointer()==array1.Get_Array_Pointer();}

    void Exchange(ARRAY_VIEW& other)
    {std::swap(m,other.m);std::swap(base_pointer,other.base_pointer);}

    void Set(T* raw_data,ID m_input)
    {m=m_input;base_pointer=raw_data;}

    void Set(const ARRAY_VIEW<T,ID>& array)
    {Set(array.base_pointer,array.m);}

    void Set(const ARRAY_VIEW<U,ID>& array)
    {Set(array.base_pointer,array.m);}

    template<class T_ARRAY>
    void Set(const ARRAY_BASE<ELEMENT,T_ARRAY,ID>& array)
    {Set(array.Derived().Get_Array_Pointer(),array.Size());}

    template<class T_ARRAY>
    void Set(ARRAY_BASE<ELEMENT,T_ARRAY,ID>& array)
    {Set(array.Derived().Get_Array_Pointer(),array.Size());}

    template<class RW>
    void Read(std::istream& input)
    {ID read_size;Read_Binary<RW>(input,read_size);
    if(read_size!=Size()){
        char buff[100];
        sprintf(buff,"Expected size %d, read %d",Value(Size()),Value(read_size));
        throw READ_ERROR(buff);}
    Read_Binary_Array<RW>(input,Get_Array_Pointer(),Value(Size()));}

    template<class RW>
    void Write(std::ostream& output) const
    {ID m=Size();Write_Binary<RW>(output,m);Write_Binary_Array<RW>(output,Get_Array_Pointer(),Value(m));}

//#####################################################################
};
template<> struct HASH_REDUCE<const char*>
{static int H(const char* key){return HASH_REDUCE<ARRAY_VIEW<const char> >::H(ARRAY_VIEW<const char>(key,(int)strlen(key)));}};
template<> struct HASH_REDUCE<char*>
{static int H(const char* key){return HASH_REDUCE<const char*>::H(key);}};
template<> struct HASH_REDUCE<std::string>
{static int H(const std::string& key){return HASH_REDUCE<ARRAY_VIEW<const char> >::H(ARRAY_VIEW<const char>(key.c_str(),(int)key.length()));}};
}
#endif
