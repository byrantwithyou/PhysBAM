//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CONVERSIONS_H
#define __CONVERSIONS_H

#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include "NUMPY.h"
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/tuple.hpp>
namespace PhysBAM{

// taken from http://language-binding.net/pyplusplus/troubleshooting_guide/automatic_conversion/automatic_conversion.html

template<class T_ARRAY>
struct ARRAY_TO_TUPLE
{
    static PyObject* convert(const T_ARRAY& array)
    {return boost::python::incref(boost::python::tuple(array).ptr());}
};

template<class T_ARRAY> struct ARRAY_IS_VECTOR:public boost::mpl::false_{};
template<class T,int d> struct ARRAY_IS_VECTOR<VECTOR<T,d> >:public boost::mpl::true_{};

template<class T_ARRAY>
struct TUPLE_TO_ARRAY
{
    typedef typename T_ARRAY::ELEMENT T;
    static const bool is_vector=ARRAY_IS_VECTOR<T_ARRAY>::value;
    static const int size_if_vector=IF<is_vector,T_ARRAY,VECTOR<int,0> >::TYPE::m;

    static bool Is_Zero(PyObject* object)
    {return PyInt_Check(object) && !PyInt_AsLong(object);}

    static void* convertible(PyObject* object)
    {if((is_vector && Is_Zero(object)) || Is_Numpy_Convertible<T_ARRAY>(object)) return object;
    if(!PySequence_Check(object)) return 0;
    try{
        boost::python::object sequence(boost::python::handle<>(boost::python::borrowed(object)));
        int size=boost::python::len(sequence);
        if(is_vector && size!=size_if_vector) return 0;
        for(int i=0;i<size;i++) if(!boost::python::extract<T>(sequence[i]).check()) return 0;
        return object;}
    catch(boost::python::error_already_set&){
        PyErr_Clear();
        return 0;}}

    static void construct_helper(PyObject* object,void* memory,boost::mpl::false_ is_vector)
    {boost::python::object sequence(boost::python::handle<>(boost::python::borrowed(object)));
    T_ARRAY* array=new (memory) T_ARRAY(boost::python::len(sequence));
    for(int i=0;i<array->Size();i++) (*array)(i+1)=boost::python::extract<T>(sequence[i])();}

    static void construct_helper(PyObject* object,void* memory,boost::mpl::true_ is_vector)
    {T_ARRAY* array=new (memory) T_ARRAY();
    if(!Is_Zero(object)){
        boost::python::object sequence(boost::python::handle<>(boost::python::borrowed(object)));
        for(int i=0;i<T_ARRAY::m;i++) (*array)[i+1]=boost::python::extract<T>(sequence[i])();}}

    static void construct(PyObject* object,boost::python::converter::rvalue_from_python_stage1_data* data)
    {typedef boost::python::converter::rvalue_from_python_storage<T_ARRAY> STORAGE;
    STORAGE* storage=reinterpret_cast<STORAGE*>(data);
    void* memory=storage->storage.bytes;
    data->convertible=memory;
    if(Is_Numpy_Array(object)){
        T_ARRAY* array=new (memory) T_ARRAY;
        Set_To_Numpy(*array,object);}
    else
        construct_helper(object,memory,boost::mpl::bool_<is_vector>());}
};

template<class T_ARRAY> void Register_Array_Conversion()
{
    boost::python::converter::registry::push_back(&TUPLE_TO_ARRAY<T_ARRAY>::convertible,&TUPLE_TO_ARRAY<T_ARRAY>::construct,boost::python::type_id<T_ARRAY>());
}

}
#endif
