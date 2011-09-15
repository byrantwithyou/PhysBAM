//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <boost/python/suite/indexing/indexing_suite.hpp>
#include <PhysBAM_Tools/Particles/PARTICLE_ATTRIBUTE.h>

namespace PhysBAM{

template<class T_ARRAY> class ARRAY_INDEXING_SUITE;
template<class T_HASH> class HASHTABLE_INDEXING_SUITE;

template<class T> struct IS_ARRAY_RESIZABLE:public TRUE_TYPE{};
template<class T,int d> struct IS_ARRAY_RESIZABLE<VECTOR<T,d> >:public FALSE_TYPE{};
template<class T> struct IS_ARRAY_RESIZABLE<RAW_ARRAY<T> >:public FALSE_TYPE{};
    
template<class T_ARRAY> 
class ARRAY_INDEXING_SUITE:public boost::python::indexing_suite<T_ARRAY,ARRAY_INDEXING_SUITE<T_ARRAY>,false,false,typename T_ARRAY::ELEMENT,int,typename T_ARRAY::ELEMENT>
{
    typedef typename T_ARRAY::ELEMENT T;
    static const bool resizable=IS_ARRAY_RESIZABLE<T_ARRAY>::value;
    typedef typename IF<resizable,T_ARRAY,ARRAY<T> >::TYPE T_SLICE_ARRAY;
    struct UNUSABLE{};

public:
    typedef T data_type;
    
    static T& get_item(T_ARRAY& array,const int i)
    {return array(i);}

    static boost::python::object get_slice(T_ARRAY& array,int start,int end)
    {start=max(start,1);end=min(end,array.Size());
    if(end<start) return boost::python::tuple();
    T_SLICE_ARRAY slice(end-start+1);for(int i=start;i<=end;i++) slice(i-start+1)=array(i);
    return boost::python::object(slice);}

    static void set_item(T_ARRAY& array,const int i,const T& v)
    {array(i)=v;}

    static void Resize_Slice(T_ARRAY& array,const int start,const typename IF<resizable,int,UNUSABLE>::TYPE end,const int new_size)
    {if(start>end+1) return;
    int old_size=end-start+1,shift=new_size-old_size;
    if(shift<0){
        for(int i=end+1;i<=array.Size();i++) array(i+shift)=array(i);
        array.Resize(array.Size()+shift);}
    else{
        int old_m=array.Size();
        array.Resize(array.Size()+shift);
        for(int i=end+1;i<=old_m;i++) array(i+shift)=array(i);}}
    
    static void Resize_Slice(T_ARRAY& array,int start,typename IF<resizable,UNUSABLE,int>::TYPE end,const int new_size)
    {if(start>end+1 || new_size==end-start+1) return;
    PyErr_SetString(PyExc_TypeError,"Can't resize a slice in a nonresizable class");
    boost::python::throw_error_already_set();}

    static void set_slice(T_ARRAY& array,int start,int end,const T& v)
    {start=max(start,1);end=min(end,array.Size());
    Resize_Slice(array,start,end,1);
    array(start)=v;}

    template<class ITERATOR>
    static void set_slice(T_ARRAY& array,int start,int end,ITERATOR first,ITERATOR last)
    {start=max(start,1);end=min(end,array.Size());
    long new_size=last-first;
    if(new_size!=(int)new_size){
        PyErr_SetString(PyExc_TypeError,"Too many elements in set_slice (int overflow)");
        boost::python::throw_error_already_set();}
    Resize_Slice(array,start,end,(int)new_size);
    for(int i=start;first!=last;++i,++first) array(i)=*first;}

    static size_t size(T_ARRAY& array)
    {return array.Size();}
    
    static int get_min_index(T_ARRAY& array)
    {return 1;}

    static int get_max_index(T_ARRAY& array)
    {return array.Size()+1;}
  
    static bool compare_index(T_ARRAY& array,const int a,const int b)
    {return a<b;}

    static int convert_index(T_ARRAY& array,PyObject* i_)
    {boost::python::extract<int> i(i_);
    if(!i.check()){
        PyErr_SetString(PyExc_TypeError, "Invalid index type");
        boost::python::throw_error_already_set();
        return 0;}
    int index=i();
    if(index<1 || index>array.Size()){
        PyErr_SetString(PyExc_IndexError, "Index out of range");
        boost::python::throw_error_already_set();}
    return index;}

    static void delete_item(T_ARRAY& array,const typename IF<resizable,int,UNUSABLE>::TYPE i)
    {for(int j=i;j<array.Size();j++) array(j)=array(j+1);
    array.Resize(array.Size()-1);}

    static void delete_item(T_ARRAY& array,const typename IF<resizable,UNUSABLE,int>::TYPE i)
    {PyErr_SetString(PyExc_TypeError,"Trying to delete an item from a nonresizable class");
    boost::python::throw_error_already_set();}

    static void delete_slice(T_ARRAY& array,int start,int end)
    {start=max(start,1);end=min(end,array.Size());
    Resize_Slice(array,start,end,0);}

    static bool contains(T_ARRAY& array,const T& key)
    {return array.Contains(key);}
};

template<class T_HASH>
class HASHTABLE_INDEXING_SUITE:public boost::python::indexing_suite<T_HASH,HASHTABLE_INDEXING_SUITE<T_HASH>,false,true,typename T_HASH::ELEMENT,typename T_HASH::KEY,typename T_HASH::KEY>
{
    typedef typename T_HASH::KEY TK;
    typedef typename T_HASH::ELEMENT T;
    struct UNUSABLE{};
public:
    typedef T data_type;
    
    static T& get_item(T_HASH& hash,TK key)
    {if(!hash.Contains(key)){ PyErr_SetString(PyExc_KeyError, "Invalid key");boost::python::throw_error_already_set();}
    return hash.Get(key);}

    static void set_item(T_HASH& hash,TK key,const T& v)
    {hash.Set(key,v);}

    // TODO: check behavior
    static void delete_item(T_HASH& hash,TK key)
    {hash.Delete(key);}
    
    static size_t size(T_HASH& hash)
    {return hash.Size();}
    
    static bool contains(T_HASH& hash,TK key)
    {return hash.Contains(key);}
  
    static TK convert_index(T_HASH& array,PyObject* i_)
    {boost::python::extract<TK> i(i_);
    if(!i.check()){
        PyErr_SetString(PyExc_TypeError, "Invalid index type");
        boost::python::throw_error_already_set();
        return 0;}
    TK index=i();
    return index;}
};
}
