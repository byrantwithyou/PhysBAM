//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCALAR_POLICY
//#####################################################################
#ifndef __SCALAR_POLICY__
#define __SCALAR_POLICY__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
#include <complex>
namespace PhysBAM{

template<class T> struct IS_SCALAR_BLOCK {static const bool value=is_scalar<T>::value;}; // true if memory layout is contiguous array of scalars
template<class T> struct IS_SCALAR_VECTOR_SPACE {static const bool value=is_scalar<T>::value;}; // true if we can compute vector space operations on the underlying array of scalars

template<class T,class ENABLER=void> struct SCALAR_POLICY{typedef struct UNUSABLE{} TYPE;};
template<class T> struct SCALAR_POLICY<T,enable_if_t<is_scalar<T>::value>>{typedef T TYPE;};
template<class T> struct SCALAR_POLICY<std::complex<T> >{typedef std::complex<T> TYPE;};
template<class T> struct SCALAR_POLICY<T,typename conditional<true,void,typename T::SCALAR>::type> {typedef typename T::SCALAR TYPE;};

}
#endif
