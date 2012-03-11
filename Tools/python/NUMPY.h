//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __NUMPY_H
#define __NUMPY_H

#define PY_ARRAY_UNIQUE_SYMBOL physbam_numpy_array_api
#ifndef PHYSBAM_IMPORT_NUMPY
#define NO_IMPORT_ARRAY
#endif

#include <boost/python/object.hpp>
#ifdef USE_NUMPY
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Point_Clouds/PARTICLES_FORWARD.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include "UTILITIES.h"
#include <numpy/arrayobject.h>
namespace PhysBAM{

//#####################################################################
// Struct NUMPY_SCALAR
//#####################################################################
// maps from primitive types to numpy type ids

struct NUMPY_NONSCALAR{};
template<class T> struct NUMPY_SCALAR:public NUMPY_NONSCALAR{};
template<> struct NUMPY_SCALAR<bool>{enum {value=PyArray_BOOL};};
template<> struct NUMPY_SCALAR<char>{enum {value=PyArray_BYTE};};
template<> struct NUMPY_SCALAR<unsigned char>{enum {value=PyArray_UBYTE};};
template<> struct NUMPY_SCALAR<short>{enum {value=PyArray_SHORT};};
template<> struct NUMPY_SCALAR<unsigned short>{enum {value=PyArray_USHORT};};
template<> struct NUMPY_SCALAR<int>{enum {value=PyArray_INT};};
template<> struct NUMPY_SCALAR<unsigned int>{enum {value=PyArray_UINT};};
template<> struct NUMPY_SCALAR<long>{enum {value=PyArray_LONG};};
template<> struct NUMPY_SCALAR<unsigned long>{enum {value=PyArray_ULONG};};
template<> struct NUMPY_SCALAR<long long>{enum {value=PyArray_LONGLONG};};
template<> struct NUMPY_SCALAR<unsigned long long>{enum {value=PyArray_ULONGLONG};};
template<> struct NUMPY_SCALAR<float>{enum {value=PyArray_FLOAT};};
template<> struct NUMPY_SCALAR<double>{enum {value=PyArray_DOUBLE};};
template<> struct NUMPY_SCALAR<long double>{enum {value=PyArray_LONGDOUBLE};};

//#####################################################################
// Struct NUMPY_INFO
//#####################################################################
// recursively extract type and shape information from statically sized types

template<class T> struct NUMPY_STATIC:public mpl::not_<boost::is_base_of<NUMPY_NONSCALAR,NUMPY_SCALAR<T> > >{};
template<class T,int d> struct NUMPY_STATIC<VECTOR<T,d> >:public mpl::true_{};
template<class T,int m,int n> struct NUMPY_STATIC<MATRIX<T,m,n> >:public mpl::true_{};

template<class T> struct NUMPY_INFO{static int Get(ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    if(!NUMPY_STATIC<T>::value){
        PyErr_SetString(PyExc_TypeError, str(boost::format("'%s' object is not convertible to numpy")%Class(T())).c_str());
        boost::python::throw_error_already_set();}
    return mpl::if_<NUMPY_STATIC<T>,NUMPY_SCALAR<T>,mpl::int_<-1> >::type::value;
}};

template<class T,int d> struct NUMPY_INFO<VECTOR<T,d> >{static int Get(ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    dimensions.Append(d);
    strides.Append(sizeof(T));
    return NUMPY_INFO<T>::Get(dimensions,strides);
}};

template<class T,int m,int n> struct NUMPY_INFO<MATRIX<T,m,n> >{static int Get(ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    dimensions.Append(m);dimensions.Append(n);
    strides.Append(sizeof(T));strides.Append(sizeof(T)*m); // note transposed ordering
    return NUMPY_INFO<T>::Get(dimensions,strides);
}};

//#####################################################################
// Function Numpy_Info
//#####################################################################
// recursively extract type and shape information from dynamically sized types

template<class TV> typename boost::enable_if<NUMPY_STATIC<TV>,int>::type
Numpy_Info(const TV& block,const void*& data,ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    data=&block;
    return NUMPY_INFO<TV>::Get(dimensions,strides);
}

template<class T_ARRAY> typename boost::enable_if<IS_ARRAY<T_ARRAY>,int>::type
Numpy_Info(const T_ARRAY& array,const void*& data,ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    typedef typename T_ARRAY::ELEMENT T;
    data=array.Get_Array_Pointer();
    dimensions.Append(array.Size());
    strides.Append(sizeof(T));
    return NUMPY_INFO<T>::Get(dimensions,strides);
}

template<class T> int
Numpy_Info(const PARTICLE_ATTRIBUTE<T>& attribute,const void*& data,ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    return Numpy_Info(attribute.array,data,dimensions,strides);
}

template<class T,class T_ARRAYS> int
Numpy_Info(const ARRAYS_ND_BASE<T,T_ARRAYS>& array,const void*& data,ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    static const npy_intp dimension=T_ARRAYS::dimension;
    const VECTOR<npy_intp,dimension> counts(array.Derived().Domain_Indices().max_corner);
    VECTOR<npy_intp,dimension> suffix_product;
    suffix_product[dimension]=sizeof(T);
    for(int i=dimension-1;i>=1;i++) suffix_product[i]=suffix_product[i+1]*counts[i+1];

    data=array.array.Get_Array_Pointer();
    dimensions.Append_Elements(counts);
    strides.Append_Elements(suffix_product);
    return NUMPY_INFO<T>::Get(dimensions,strides);
}

template<class T> int
Numpy_Info(const MATRIX_MXN<T>& matrix,const void*& data,ARRAY<npy_intp>& dimensions,ARRAY<npy_intp>& strides)
{
    data=matrix.x;
    dimensions.Append(matrix.Rows());dimensions.Append(matrix.Columns());
    strides.Append(sizeof(T));strides.Append(sizeof(T)*matrix.Rows()); // note transposed ordering
    return NUMPY_INFO<T>::Get(dimensions,strides);
}

//#####################################################################
// Function Numpy_Shape_Match
//#####################################################################
// check whether dynamic type can be resized to fit a given numpy array

template<class TV> typename boost::enable_if<NUMPY_STATIC<TV>,bool>::type
Numpy_Shape_Match(FIRST<TV>,RAW_ARRAY<const npy_intp> dimensions)
{
    ARRAY<npy_intp> subdimensions,strides;
    NUMPY_INFO<TV>::Get(subdimensions,strides);
    return dimensions==subdimensions;
}

template<class T_ARRAY> typename boost::enable_if<IS_ARRAY<T_ARRAY>,bool>::type
Numpy_Shape_Match(FIRST<T_ARRAY>,RAW_ARRAY<const npy_intp> dimensions)
{
    typedef typename T_ARRAY::ELEMENT T;
    ARRAY<npy_intp> subdimensions,strides;
    NUMPY_INFO<T>::Get(subdimensions,strides);
    int rank=dimensions.Size(),subrank=subdimensions.Size();
    return rank==subrank+1 && dimensions.Subset(IDENTITY_MAP<>(subrank)+1)==subdimensions;
}

template<class T,class T_ARRAYS> bool
Numpy_Shape_Match(FIRST<ARRAYS_ND_BASE<T,T_ARRAYS> >,RAW_ARRAY<const npy_intp> dimensions)
{
    static const int dimension=T_ARRAYS::dimension;
    ARRAY<npy_intp> subdimensions,strides;
    NUMPY_INFO<T>::Get(subdimensions,strides);
    int rank=dimensions.Size(),subrank=subdimensions.Size();
    return rank==subrank+dimension && dimensions.Subset(IDENTITY_MAP<>(subrank)+dimension)==subdimensions;
}
template<class T> bool Numpy_Shape_Match(FIRST<ARRAYS<VECTOR<T,1> > >,RAW_ARRAY<const npy_intp> dimensions){return Numpy_Shape_Match(FIRST<ARRAYS_ND_BASE<T,ARRAYS<VECTOR<T,1> > > >(),dimensions);}
template<class T> bool Numpy_Shape_Match(FIRST<ARRAYS<VECTOR<T,2> > >,RAW_ARRAY<const npy_intp> dimensions){return Numpy_Shape_Match(FIRST<ARRAYS_ND_BASE<T,ARRAYS<VECTOR<T,2> > > >(),dimensions);}
template<class T> bool Numpy_Shape_Match(FIRST<ARRAYS<VECTOR<T,3> > >,RAW_ARRAY<const npy_intp> dimensions){return Numpy_Shape_Match(FIRST<ARRAYS_ND_BASE<T,ARRAYS<VECTOR<T,3> > > >(),dimensions);}

template<class T> bool
Numpy_Shape_Match(FIRST<MATRIX_MXN<T> >,RAW_ARRAY<const npy_intp> dimensions)
{
    ARRAY<npy_intp> subdimensions,strides;
    NUMPY_INFO<T>::Get(subdimensions,strides);
    int rank=dimensions.Size(),subrank=subdimensions.Size();
    return rank==subrank+2 && dimensions.Subset(IDENTITY_MAP<>(subrank)+2)==subdimensions;
}

//#####################################################################
// Function Numpy_Resize
//#####################################################################
// resize a physbam array to match a given numpy array

template<class TV> typename boost::enable_if<NUMPY_STATIC<TV> >::type
Numpy_Resize(TV& block,RAW_ARRAY<const npy_intp> dimensions)
{
    // static types needn't be resized
}

template<class T_ARRAY> typename boost::enable_if<IS_ARRAY<T_ARRAY> >::type
Numpy_Resize(T_ARRAY& array,RAW_ARRAY<const npy_intp> dimensions)
{
    PHYSBAM_ASSERT(dimensions.Size()>=1);
    array.Resize(dimensions(1),false,false);
}

template<class T,class T_ARRAYS> void
Numpy_Resize(const ARRAYS_ND_BASE<T,T_ARRAYS>& array,RAW_ARRAY<const npy_intp> dimensions)
{
    static const int dimension=T_ARRAYS::dimension;
    PHYSBAM_ASSERT(dimensions.Size()>=dimension);
    BOX<VECTOR<int,dimension> > indices;
    for(int a=0;a<dimension;a++){
        indices.min_corner[a]=1;
        indices.max_corner[a]=dimensions(a);}
    array.Resize(indices,false,false);
}

template<class T> void
Numpy_Resize(MATRIX_MXN<T>& matrix,RAW_ARRAY<const npy_intp> dimensions)
{
    PHYSBAM_ASSERT(dimensions.Size()>=2);
    matrix.Resize(dimensions(1),dimensions(2));
}

//#####################################################################
// Function Numpy_View
//#####################################################################
template<class T_ARRAY> boost::python::handle<> Numpy_View(T_ARRAY& array)
{
    // extract memory layout information
    const void* data;
    ARRAY<npy_intp> dimensions,strides;
    int numpy_type=Numpy_Info(array,data,dimensions,strides);
    PHYSBAM_ASSERT(dimensions.Size()==strides.Size());

    // wrap the existing array as a numpy array without copying data
    int flags=0;
    if(!IS_CONST<T_ARRAY>::value) flags|=NPY_WRITEABLE;
    return boost::python::handle<>(PyArray_NewFromDescr(&PyArray_Type,PyArray_DescrFromType(numpy_type),
        dimensions.Size(),dimensions.Get_Array_Pointer(),strides.Get_Array_Pointer(),
        const_cast<void*>(data),flags,0));
}
//#####################################################################
// Function As_Numpy
//#####################################################################
template<class T_ARRAY> boost::python::handle<> As_Numpy(const T_ARRAY& array)
{
    // copy the numpy array before anyone has a chance to resize it
    return boost::python::handle<>(PyArray_FromAny(Numpy_View(array).get(),0,0,INT_MAX,NPY_ENSURECOPY,Py_None));
}

template<class T_ARRAY> boost::python::handle<> As_Numpy_With_Context(const T_ARRAY& array,const boost::python::object& context)
{
    return As_Numpy(array); // ignore context
}
//#####################################################################
// Function Numpy_Dimensions
//#####################################################################
inline RAW_ARRAY<const npy_intp> Numpy_Dimensions(PyObject* object)
{
    return RAW_ARRAY<const npy_intp>(PyArray_NDIM(object),PyArray_DIMS(object));
}
//#####################################################################
// Function Is_Numpy_Array
//#####################################################################
inline bool Is_Numpy_Array(PyObject* object)
{
    return PyArray_Check(object);
}
//#####################################################################
// Function Is_Numpy_Convertible
//#####################################################################
template<class T_ARRAY> bool Is_Numpy_Convertible(PyObject* object)
{
    return PyArray_Check(object) && Numpy_Shape_Match(FIRST<T_ARRAY>(),Numpy_Dimensions(object));
}
//#####################################################################
// Function Set_To_Numpy
//#####################################################################
template<class T_ARRAY> void Set_To_Numpy(T_ARRAY& array,PyObject* object)
{
    PHYSBAM_ASSERT(PyArray_Check(object));
    RAW_ARRAY<const npy_intp> dimensions(Numpy_Dimensions(object));

    // resize array
    Numpy_Resize(array,dimensions);

    // fill it with data
    boost::python::handle<> wrapper=Numpy_View(array);
    if(PyArray_CastTo((PyArrayObject*)wrapper.get(),(PyArrayObject*)object)<0){
        boost::python::object from(boost::python::handle<>(boost::python::borrowed(object))),to(wrapper);
        PyErr_SetString(PyExc_ValueError,str(boost::format("Conversion from %s (shape = %s, dtype = %s) to %s (shape = %s, dtype = %s) failed")
            %Class(from)%Repr(from.attr("shape"))%Repr(from.attr("dtype").attr("name"))
            %Class(array)%Repr(to.attr("shape"))%Repr(to.attr("dtype").attr("name"))).c_str());
        boost::python::throw_error_already_set();}
}
//#####################################################################
// Function Import_Numpy
//#####################################################################
#ifdef PHYSBAM_IMPORT_NUMPY
void Import_Numpy()
{
    if(_import_array()<0){
        PyErr_Print();
        PyErr_SetString(PyExc_ImportError,"numpy.core.multiarray failed to import");
        boost::python::throw_error_already_set();}
}
#endif
//#####################################################################
}
#else
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T_ARRAY> boost::python::object As_Numpy(const T_ARRAY& array)
{
    PHYSBAM_NOT_IMPLEMENTED("Conversion to numpy");
}

template<class T_ARRAY> boost::python::object As_Numpy_With_Context(const T_ARRAY& array,const boost::python::object& context)
{
    PHYSBAM_NOT_IMPLEMENTED("Conversion to numpy");
}

inline bool Is_Numpy_Array(PyObject* object)
{
    return false; // can't detect numpy with numpy library
}

template<class T_ARRAY> bool Is_Numpy_Convertible(PyObject* object)
{
    return false; // can't detect numpy with numpy library
}

template<class T_ARRAY> void Set_To_Numpy(T_ARRAY& array,PyObject* object)
{
    PHYSBAM_NOT_IMPLEMENTED("Conversion from numpy");
}

#ifdef PHYSBAM_IMPORT_NUMPY
void Import_Numpy()
{}
#endif

}
#endif
#endif
