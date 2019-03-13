//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY
//#####################################################################
#ifndef __ARRAY__
#define __ARRAY__

#include <Core/Arrays/ARRAY_BASE.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/exchange.h>
#include <Core/Math_Tools/min.h>
#include <Core/Utilities/EXCEPTIONS.h>
#include <Core/Utilities/PHYSBAM_ATTRIBUTE.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <initializer_list>
namespace PhysBAM{

template<class T,class ID> class ARRAY;
template<class T,class ID> struct IS_ARRAY<ARRAY<T,ID> > {static const bool value=true;};
template<class T,class ID> struct IS_ARRAY_VIEW<ARRAY<T,ID> > {static const bool value=false;};

template<class T,class ID> struct CANONICALIZE_CONST_ARRAY<ARRAY<T,ID> >:public FIRST<ARRAY_VIEW<const T,ID> >{};

struct NO_INIT{};
struct USE_INIT{};
struct INIT_ALL{};
static const NO_INIT no_init={};
static const USE_INIT use_init={};
static const INIT_ALL init_all={};

template<class T,class ID>
class ARRAY:public ARRAY_BASE<T,ARRAY<T,ID>,ID>
{
    struct UNUSABLE{};
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    template<class T2> struct REBIND{typedef ARRAY<T2,ID> TYPE;};
    typedef T ELEMENT;typedef ID INDEX;
    using ARRAY_BASE<T,ARRAY<T,ID>,ID>::Same_Array;
    typedef T* iterator; // for stl
    typedef const T* const_iterator; // for stl
    T* base_pointer;
    ID buffer_size;

    ID m; // the current size of the array (buffer_size may be larger for elbow room)

    ARRAY()
        :base_pointer(0),buffer_size(0),m(0)
    {}

    explicit ARRAY(ID m_input)
        :base_pointer(0),buffer_size(m_input),m(m_input)
    {
        assert(m>=ID());
        base_pointer=(T*)new unsigned char[sizeof(T)*Value(m)];
        for(int i=0,j=Value(m);i<j;i++) new(base_pointer+i)T();
    }

    ARRAY(ID m_input,NO_INIT)
        :base_pointer(0),buffer_size(m_input),m(m_input)
    {
        assert(m>=ID());
        base_pointer=(T*)new unsigned char[sizeof(T)*Value(m)];
        for(int i=0,j=Value(m);i<j;i++) new(base_pointer+i)T;
    }

    ARRAY(ID m_input,USE_INIT,const T& __restrict__ initialization_value)
        :buffer_size(m_input),m(m_input)
    {
        assert(m>=ID());
        base_pointer=(T*)new unsigned char[sizeof(T)*Value(m)];
        for(int i=0,j=Value(m);i<j;i++) new(base_pointer+i)T(initialization_value);
    }

    explicit ARRAY(const INITIAL_SIZE m_input)
        :ARRAY(ID(Value(m_input)))
    {
    }

    ARRAY(const INITIAL_SIZE m_input,NO_INIT)
        :ARRAY(ID(Value(m_input)),no_init)
    {
    }

    ARRAY(const INITIAL_SIZE m_input,USE_INIT,const T& __restrict__ initialization_value)
        :ARRAY(ID(Value(m_input)),use_init,initialization_value)
    {
    }

    ARRAY(const ARRAY& array)
        :buffer_size(array.m),m(array.m)
    {
        base_pointer=(T*)new unsigned char[sizeof(T)*Value(array.m)];
        for(int i=0,j=Value(m);i<j;i++) new(base_pointer+i)T(array.base_pointer[i]);
    }

    ARRAY(ARRAY&& array)
        :base_pointer(array.base_pointer),buffer_size(array.buffer_size),m(array.m)
    {
        array.base_pointer=0;
        array.buffer_size=ID();
        array.m=ID();
    }

    template<class T_ARRAY>
    explicit ARRAY(const ARRAY_BASE<T,T_ARRAY,ID>& array)
        :buffer_size(array.Size()),m(array.Size())
    {
        base_pointer=(T*)new unsigned char[sizeof(T)*Value(m)];
        for(int i=0,j=Value(m);i<j;i++) new(base_pointer+i)T(array(ID(i)));
    }

    ARRAY(std::initializer_list<T> init)
        :buffer_size(init.size()),m(init.size())
    {
        base_pointer=(T*)new unsigned char[sizeof(T)*Value(m)];
        T* p=base_pointer;
        for(const auto& a:init) new(p++)T(a);
    }

    ~ARRAY()
    {
        Call_Destructors_And_Free();
    }

    ID Size() const
    {return m;}

    T& operator()(ID i)
    {assert((unsigned)Value(i)<(unsigned)Value(m));return base_pointer[Value(i)];}

    const T& operator()(ID i) const
    {assert((unsigned)Value(i)<(unsigned)Value(m));return base_pointer[Value(i)];}

    ARRAY& operator=(const ARRAY& array)
    {
        if(buffer_size<array.m){
            Call_Destructors_And_Free();
            buffer_size=m=array.m;
            base_pointer=(T*)new unsigned char[sizeof(T)*Value(array.m)];
            for(int i=0,j=Value(m);i<j;i++) new(base_pointer+i)T(array.base_pointer[i]);
            return *this;}
        if(Same_Array(*this,array)) return *this;
        int n=Value(m),p=Value(array.m);
        for(int i=0,j=std::min(n,p);i<j;i++) base_pointer[i]=array.base_pointer[i];
        for(int i=n;i<p;i++) new(base_pointer+i)T(array.base_pointer[i]);
        for(int i=p;i<n;i++) base_pointer[i].~T();
        m=array.m;
        return *this;
    }

    ARRAY& operator=(ARRAY&& array)
    {
        Call_Destructors_And_Free();
        base_pointer=array.base_pointer;
        buffer_size=array.buffer_size;
        m=array.m;
        array.base_pointer=0;
        array.buffer_size=ID();
        array.m=ID();
        return *this;
    }

    template<class T_ARRAY>
    ARRAY& operator=(const ARRAY_BASE<T,T_ARRAY,ID>& array)
    {
        ID array_m=array.Size();
        if(buffer_size<array_m){
            Call_Destructors_And_Free();
            buffer_size=m=array.Size();
            base_pointer=(T*)new unsigned char[sizeof(T)*Value(m)];
            for(int i=0,j=Value(m);i<j;i++) new(base_pointer+i)T(array(ID(i)));
            return *this;}
        int n=Value(m),p=Value(array_m);
        for(int i=0,j=std::min(n,p);i<j;i++) base_pointer[i]=array(ID(i));
        for(int i=n;i<p;i++) new(base_pointer+i)T(array(ID(i)));
        for(int i=p;i<n;i++) base_pointer[i].~T();
        m=array_m;
        return *this;
    }

    bool Valid_Index(ID i) const
    {return (unsigned)Value(i)<(unsigned)Value(m);}

    T* Get_Array_Pointer()
    {return base_pointer;}

    const T* Get_Array_Pointer() const
    {return base_pointer;}

    ID Max_Size() const
    {return buffer_size;}

    void Compact()
    {Exact_Resize(m);}

private:
    void Call_Destructors_And_Free() PHYSBAM_ALWAYS_INLINE
    {
        for(int i=Value(m)-1;i>=0;i--) base_pointer[i].~T();
        delete[] (unsigned char*)base_pointer;
    }

    void Resize_Helper_Fill(int n) PHYSBAM_ALWAYS_INLINE
    {
        for(int i=0;i<n;i++) base_pointer[i]=T();
    }

    void Resize_Helper_Fill(int n,const T& initialization_value) PHYSBAM_ALWAYS_INLINE
    {
        for(int i=0;i<n;i++) base_pointer[i]=initialization_value;
    }

    template<bool exact,bool copy,bool init,bool fill,class... Args>
    void Resize_Helper(ID m_new,Args&&... initialization_value)
    {
        int n=Value(m),p=Value(m_new);
        if(buffer_size<m_new){
            int new_buffer_size=exact?p:(p?2*p:2);
            T* new_buffer=(T*)new unsigned char[sizeof(T)*new_buffer_size];
            if(init) for(int i=(fill?0:n);i<p;i++) new(new_buffer+i)T(initialization_value...);
            if(copy) for(int i=0;i<n;i++) new(new_buffer+i)T(std::move(base_pointer[i]));
            if(!init) for(int i=(copy?n:0);i<p;i++) new(new_buffer+i)T;
            Call_Destructors_And_Free();
            base_pointer=new_buffer;
            buffer_size=ID(new_buffer_size);}
        else{
            if(fill) Resize_Helper_Fill(std::min(n,p),initialization_value...);
            if(init) for(int i=n;i<p;i++) new(base_pointer+i)T(initialization_value...);
            else for(int i=n;i<p;i++) new(base_pointer+i)T;
            for(int i=p;i<n;i++) base_pointer[i].~T();}
        m=m_new;
    }

    template<class... Args>
    void Append_Helper(Args&&... initialization_value)
    {
        int n=Value(m);
        int new_buffer_size=(n?2*n:2);
        T* new_buffer=(T*)new unsigned char[sizeof(T)*new_buffer_size];
        for(int i=0;i<n;i++) new(new_buffer+i)T(std::move(base_pointer[i]));
        new(new_buffer+n)T(std::forward<Args>(initialization_value)...);
        Call_Destructors_And_Free();
        base_pointer=new_buffer;
        buffer_size=ID(new_buffer_size);
    }

public:
    void Preallocate(ID max_size)
    {
        if(buffer_size<max_size){
            int n=Value(m);
            T* new_buffer=(T*)new unsigned char[sizeof(T)*Value(max_size)];
            for(int i=0;i<n;i++) new(new_buffer+i)T(std::move(base_pointer[i]));
            Call_Destructors_And_Free();
            base_pointer=new_buffer;
            buffer_size=max_size;}
    }

    void Resize(ID m_new)
    {Resize_Helper<false,true,true,false>(m_new);}
    
    void Resize(ID m_new,NO_INIT)
    {Resize_Helper<false,false,false,false>(m_new);}

    // Aliasing with initialization_value is okay.
    void Resize(ID m_new,USE_INIT,const T& initialization_value)
    {Resize_Helper<false,true,true,false>(m_new,initialization_value);}

    void Resize(ID m_new,INIT_ALL)
    {Resize_Helper<false,false,true,true>(m_new);}

    // Aliasing with initialization_value is okay.
    void Resize(ID m_new,INIT_ALL,const T& initialization_value)
    {Resize_Helper<false,false,true,true>(m_new,initialization_value);}

    void Exact_Resize(ID m_new)
    {Resize_Helper<true,true,true,false>(m_new);}
    
    void Exact_Resize(ID m_new,NO_INIT)
    {Resize_Helper<true,false,false,false>(m_new);}

    // Aliasing with initialization_value is okay.
    void Exact_Resize(ID m_new,USE_INIT,const T& initialization_value)
    {Resize_Helper<true,true,true,false>(m_new,initialization_value);}

    void Exact_Resize(ID m_new,INIT_ALL)
    {Resize_Helper<true,false,true,true>(m_new);}

    // Aliasing with initialization_value is okay.
    void Exact_Resize(ID m_new,INIT_ALL,const T& initialization_value)
    {Resize_Helper<true,false,true,true>(m_new,initialization_value);}

    // Aliasing with element is okay.
    ID Append(const T& element) PHYSBAM_ALWAYS_INLINE
    {if(buffer_size<=m) Append_Helper(T(element));else new(base_pointer+Value(m))T(element);return m++;}

    // Aliasing with element is okay.
    ID Append(T&& element) PHYSBAM_ALWAYS_INLINE
    {if(buffer_size<=m) Append_Helper(std::move(element));else new(base_pointer+Value(m))T(std::move(element));return m++;}

    // Aliasing is not allowed
    template<class T_ARRAY>
    void Append_Elements(const ARRAY_BASE<T,T_ARRAY,ID>& append_array)
    {ID n=append_array.Size(),p=Value(m);Resize_Helper<false,true,false,false>(m+n);for(ID i(0);i<n;i++) (*this)(i+p)=append_array(i);}

    // Note: not very efficient
    void Append_Unique(const T& element)
    {for(ID i(0);i<m;i++) if((*this)(i)==element) return;Append(element);}

    // Note: not efficient
    template<class T_ARRAY>
    void Append_Unique_Elements(const ARRAY_BASE<T,T_ARRAY,ID>& append_array)
    {ID n=append_array.Size();for(ID i(0);i<n;i++) Append_Unique(append_array(i));}

    ID Add_End() PHYSBAM_ALWAYS_INLINE
    {if(buffer_size<=m) Append_Helper();else new(base_pointer+Value(m))T();return m++;}

    void Remove_End()
    {assert(m>ID());base_pointer[Value(--m)].~T();}

    void Remove_Index(ID index) // preserves ordering of remaining elements
    {assert((unsigned)Value(index)<(unsigned)Value(m));for(ID i=index;i<m-1;i++) (*this)(i)=std::move((*this)(i+1));Remove_End();}

    void Remove_Index_Lazy(ID index)
    {assert((unsigned)Value(index)<(unsigned)Value(m));if(index<m-1) (*this)(index)=std::move((*this)(m-1));Remove_End();}

    void Remove_All() // if elements are non-primitive this may waste memory
    {for(int i=Value(m)-1;i>=0;i--) base_pointer[i].~T();m=ID();}

    void Clean_Memory()
    {Call_Destructors_And_Free();base_pointer=0;buffer_size=m=ID();}

    void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
    {for(ID i(0);i<m;i++) delete (*this)(i);Clean_Memory();}

    // Aliasing is not allowed
    void Insert(const T& element,ID index)
    {Add_End();for(ID i=m-1;i>index;i--) (*this)(i)=std::move((*this)(i-1));(*this)(index)=element;}

    T Pop_Value()
    {T copy=base_pointer[Value(m)-1];Remove_End();return copy;}

    void Pop()
    {Remove_End();}

    void Pop_Elements(const int count) // return value should be copied immediately, not kept around
    {int n=Value(m);assert(n-count>=0);for(int i=n-1;i>=n-count;i--) base_pointer[i].~T();m-=count;}

    T* begin() // for stl
    {return Get_Array_Pointer();}

    const T* begin() const // for stl
    {return Get_Array_Pointer();}

    T* end() // for stl
    {return Get_Array_Pointer()+Value(m);}

    const T* end() const // for stl
    {return Get_Array_Pointer()+Value(m);}

    void Exchange(ARRAY<T,ID>& other)
    {exchange(m,other.m);exchange(base_pointer,other.base_pointer);exchange(buffer_size,other.buffer_size);}

    template<class RW> void Read(std::istream& input)
    {Clean_Memory();ID m;Read_Binary<RW>(input,m);
    if(m<ID()){
        char buff[100];
        sprintf(buff,"Invalid negative array size %d",Value(m));
        throw READ_ERROR(buff);}
    if(!m) return;
    Exact_Resize(m);
    Read_Binary_Array<RW>(input,Get_Array_Pointer(),Value(m));}

    template<class RW> void Write(std::ostream& output) const
    {Write_Prefix<RW>(output,m);}

    template<class RW> void Write_Prefix(std::ostream& output,ID prefix) const
    {PHYSBAM_ASSERT(ID(0)<=prefix && prefix<=Size());
    Write_Binary<RW>(output,Value(prefix));Write_Binary_Array<RW>(output,Get_Array_Pointer(),Value(prefix));}
//#####################################################################
};
}
#endif
