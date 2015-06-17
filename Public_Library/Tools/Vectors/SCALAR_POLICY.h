//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCALAR_POLICY
//#####################################################################
#ifndef __SCALAR_POLICY__
#define __SCALAR_POLICY__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T> struct is_scalar_BLOCK {static const bool value=is_scalar<T>::value;}; // true if memory layout is contiguous array of scalars
template<class T> struct is_scalar_VECTOR_SPACE {static const bool value=is_scalar<T>::value;}; // true if we can compute vector space operations on the underlying array of scalars

template<class T,class ENABLER=void> struct SCALAR_POLICY{typedef struct UNUSABLE{} TYPE;};
template<class T> struct SCALAR_POLICY<T,typename ENABLE_IF<is_scalar<T>::value>::TYPE>{typedef T TYPE;};
template<class T> struct SCALAR_POLICY<T,typename IF<true,void,typename T::SCALAR>::TYPE> {typedef typename T::SCALAR TYPE;};

}
#endif
