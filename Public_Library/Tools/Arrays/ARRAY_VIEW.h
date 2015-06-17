//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_VIEW
//#####################################################################
#ifndef __ARRAY_VIEW__
#define __ARRAY_VIEW__

#include <Tools/Arrays/ARRAY_BASE.h>
#include <Tools/Log/LOG_PRINTF.h>
#include <Tools/Utilities/EXCEPTIONS.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class ID> struct IS_ARRAY<ARRAY_VIEW<T,ID> > {static const bool value=true;};
template<class T,class ID> struct IS_ARRAY_VIEW<ARRAY_VIEW<T,ID> > {static const bool value=true;};
template<class T,class ID,class T_NEW> struct REBIND<ARRAY_VIEW<T,ID>,T_NEW>{typedef ARRAY_VIEW<typename conditional<is_const<T>::value,const T_NEW,T_NEW>::type,ID> TYPE;};

template<class T,class ID>
class ARRAY_VIEW:public ARRAY_BASE<typename remove_const<T>::type,ARRAY_VIEW<T,ID>,ID>
{
    struct UNUSABLE{};
    template<class S> struct COPY_CONST{typedef typename conditional<is_const<T>::value,typename add_const<S>::type,S>::type TYPE;};
    typedef ARRAY_BASE<typename remove_const<T>::type,ARRAY_VIEW<T,ID>,ID> BASE;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef typename remove_const<T>::type ELEMENT;typedef ID INDEX;
    typedef T& RESULT_TYPE;

    // m and base_pointer inherit constness of T
    typename COPY_CONST<ID>::TYPE m;
    friend class ARRAY_VIEW<typename conditional<is_const<T>::value,ELEMENT,const ELEMENT>::type,ID>;
    typename COPY_CONST<T*>::TYPE base_pointer;

    using BASE::Same_Array;

    ARRAY_VIEW()
        :m(0),base_pointer(0)
    {}

    ARRAY_VIEW(const ID m,T* raw_data)
        :m(m),base_pointer(raw_data)
    {}

    ARRAY_VIEW(const ARRAY_VIEW<T,ID>& array)
        :m(array.m),base_pointer(array.base_pointer)
    {}

    ARRAY_VIEW(const ARRAY_VIEW<typename conditional<is_const<T>::value,typename remove_const<T>::type,UNUSABLE>::type,ID>& array)
        :m(array.m),base_pointer(array.base_pointer)
    {}

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY& array,typename enable_if<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && !IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::type unusable=UNUSABLE())
        :m(array.Size()),base_pointer(array.Get_Array_Pointer())
    {}

    template<class T_ARRAY>
    ARRAY_VIEW(T_ARRAY array,typename enable_if<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::type unusable=UNUSABLE())
        :m(array.Size()),base_pointer(array.Get_Array_Pointer())
    {}

    ID Size() const
    {return m;}

    T& operator()(const ID i)
    {assert((unsigned)Value(i)<(unsigned)Value(m));return base_pointer[Value(i)];}

    const T& operator()(const ID i) const
    {assert((unsigned)Value(i)<(unsigned)Value(m));return base_pointer[Value(i)];}

    bool Valid_Index(const ID i) const
    {return (unsigned)Value(i)<(unsigned)Value(m);}

    ARRAY_VIEW& operator=(const ARRAY_VIEW& source)
    {return BASE::operator=(source);}

    template<class T_ARRAY1>
    ARRAY_VIEW& operator=(const T_ARRAY1& source)
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

    ARRAY_VIEW<typename remove_const<T>::type>& Const_Cast() const // return reference to allow Exchange
    {return reinterpret_cast<ARRAY_VIEW<typename remove_const<T>::type>&>(const_cast<ARRAY_VIEW&>(*this));}

    static bool Same_Array(const ARRAY_VIEW& array0,const ARRAY_VIEW& array1)
    {return array0.Get_Array_Pointer()==array1.Get_Array_Pointer();}

private:
    template<class T2>
    static void exchange(T2& a,T2& b)
    {T2 tmp=a;a=b;b=tmp;}

public:
    void Exchange(ARRAY_VIEW& other)
    {STATIC_ASSERT(!is_const<T>::value); // make ARRAY_VIEW<const T> equivalent to const ARRAY_VIEW<const T>
    exchange(m,other.m);exchange(base_pointer,other.base_pointer);}

    void Set(const ID m_input,typename COPY_CONST<T*>::TYPE raw_data)
    {m=m_input;base_pointer=raw_data;}

    void Set(const ARRAY_VIEW<T,ID>& array)
    {Set(array.m,array.base_pointer);}

    void Set(const ARRAY_VIEW<typename conditional<is_const<T>::value,typename remove_const<T>::type,UNUSABLE>::type,ID>& array)
    {Set(array.m,array.base_pointer);}

    template<class T_ARRAY>
    void Set(T_ARRAY& array,typename enable_if<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && !IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::type unusable=UNUSABLE())
    {Set(array.m,array.Get_Array_Pointer());}

    template<class T_ARRAY>
    void Set(T_ARRAY array,typename enable_if<is_same<ELEMENT,typename T_ARRAY::ELEMENT>::value && IS_ARRAY_VIEW<T_ARRAY>::value,UNUSABLE>::type unusable=UNUSABLE())
    {Set(array.m,array.Get_Array_Pointer());}

    template<class RW>
    void Read(std::istream& input)
    {ID read_size;Read_Binary<RW>(input,read_size);
    if(read_size!=Size()) throw READ_ERROR(LOG::sprintf("Expected size %d, read %d",Size(),read_size));
    Read_Binary_Array<RW>(input,Get_Array_Pointer(),Value(Size()));}

    template<class RW>
    void Write(std::ostream& output) const
    {ID m=Size();Write_Binary<RW>(output,m);Write_Binary_Array<RW>(output,Get_Array_Pointer(),Value(m));}

//#####################################################################
};
template<> struct HASH_REDUCE<const char*>
{static int H(const char* key){return HASH_REDUCE<ARRAY_VIEW<const char> >::H(ARRAY_VIEW<const char>((int)strlen(key),key));}};
template<> struct HASH_REDUCE<char*>
{static int H(const char* key){return HASH_REDUCE<const char*>::H(key);}};
template<> struct HASH_REDUCE<std::string>
{static int H(const std::string& key){return HASH_REDUCE<ARRAY_VIEW<const char> >::H(ARRAY_VIEW<const char>((int)key.length(),key.c_str()));}};
}
#endif
