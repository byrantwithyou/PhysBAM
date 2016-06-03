//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ARRAYS_FORWARD
//#####################################################################
#ifndef __ARRAYS_FORWARD__
#define __ARRAYS_FORWARD__

#include <Tools/Utilities/STATIC_ASSERT.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,class ID=int> class ARRAY_VIEW;
template<class T,class T_ARRAY,class ID=int> class ARRAY_BASE;
template<class T,class T_ARRAY,class ID=int> class ARRAY_EXPRESSION;
template<class T,class T_ARRAY,class ID=int> class INDEX_ITERATED_ARRAY;

template<class T,class ID=int> class ARRAY;

template<class ID=int> class IDENTITY_ARRAY;
template<class T,class ID=int> class CONSTANT_ARRAY;

template<class T_ARRAY,class T_INDICES=ARRAY<int>&> class INDIRECT_ARRAY;

template<class T_ARRAY,class T_PROJECTOR> class PROJECTED_ARRAY;
template<class T_STRUCT,class T_FIELD,T_FIELD T_STRUCT::* field> struct FIELD_PROJECTOR;
struct INDEX_PROJECTOR;

template<class T,int d> class VECTOR;

template<class T_ARRAY> struct ARRAY_RESULT_TYPE{typedef typename T_ARRAY::RESULT_TYPE TYPE;};
template<class T_ARRAY> struct ARRAY_RESULT_TYPE<const T_ARRAY>{typedef typename T_ARRAY::CONST_RESULT_TYPE TYPE;};

template<class T,class ID,class SCALAR> struct REPLACE_FLOATING_POINT<ARRAY<T,ID>,SCALAR>{typedef ARRAY<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,ID> TYPE;};
template<class T_ARRAY0,class T_ARRAY1> class ARRAY_SUM;
template<class T_ARRAY0,class T_ARRAY1> class ARRAY_DIFFERENCE;
template<class T_ARRAY0,class T_ARRAY1> class ARRAY_PRODUCT;
template<class T_ARRAY0,class T_ARRAY1> class ARRAY_QUOTIENT;
template<class T_ARRAY> class ARRAY_NEGATION;
template<class T0,class T_ARRAY1> class ARRAY_PLUS_SCALAR;
template<class T0,class T_ARRAY1> class ARRAY_LEFT_MULTIPLE;

template<class TV> class RANGE;
template<int d> class FACE_INDEX;
template<class T> struct DOMAIN_INDEX_TYPE {typedef INTERVAL<T> TYPE;};
template<int d> struct DOMAIN_INDEX_TYPE<VECTOR<int,d> > {typedef RANGE<VECTOR<int,d> > TYPE;};
template<int d> struct DOMAIN_INDEX_TYPE<FACE_INDEX<d> > {typedef RANGE<VECTOR<int,d> > TYPE;};


}
#endif
