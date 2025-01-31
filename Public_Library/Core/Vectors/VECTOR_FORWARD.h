//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header VECTOR_FORWARD
//#####################################################################
#ifndef __VECTOR_FORWARD__
#define __VECTOR_FORWARD__

#include <Core/Read_Write/READ_WRITE_FORWARD.h>
#include <Core/Utilities/STATIC_ASSERT.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;
template<class T,int d> class ZERO_VECTOR;

template<class T> class INTERVAL;
template<class T> class QUATERNION;
template<class TV> class ROTATION;
template<class T> class TETRAHEDRAL_GROUP;
template<class TV> class FRAME;
template<class TV> class TWIST;

template<class T_VECTOR1,class T_VECTOR2> class VECTOR_SUM;
template<class T_VECTOR1,class T_VECTOR2> class VECTOR_DIFFERENCE;
template<class T_VECTOR,class T2> class VECTOR_SCALE;
template<class T_VECTOR> class VECTOR_NEGATION;

template<class T> struct IS_SCALAR_BLOCK;
template<class T> struct IS_SCALAR_VECTOR_SPACE;
template<class T,int d> struct IS_SCALAR_BLOCK<VECTOR<T,d> > {static const bool value=(d>0) && IS_SCALAR_BLOCK<T>::value;};
template<class T,int d> struct IS_SCALAR_VECTOR_SPACE<VECTOR<T,d> > {static const bool value=(d>0) && IS_SCALAR_VECTOR_SPACE<T>::value;};
template<class T,int d,class RW> struct IS_BINARY_IO_SAFE<VECTOR<T,d>,RW> {static const bool value=(d>0) && IS_BINARY_IO_SAFE<T,RW>::value;};

template<class T> struct HAS_CHEAP_COPY;
template<class T,int d> struct HAS_CHEAP_COPY<VECTOR<T,d> > {static const bool value=true;};

template<class T,int d,class SCALAR> struct REPLACE_FLOATING_POINT<VECTOR<T,d>,SCALAR>{typedef VECTOR<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,d> TYPE;};
template<class TV,class SCALAR> struct REPLACE_FLOATING_POINT<TWIST<TV>,SCALAR>{typedef TWIST<typename REPLACE_FLOATING_POINT<TV,SCALAR>::TYPE> TYPE;};
template<class TV,class SCALAR> struct REPLACE_FLOATING_POINT<FRAME<TV>,SCALAR>{typedef FRAME<typename REPLACE_FLOATING_POINT<TV,SCALAR>::TYPE> TYPE;};
template<class TV,class SCALAR> struct REPLACE_FLOATING_POINT<ROTATION<TV>,SCALAR>{typedef ROTATION<typename REPLACE_FLOATING_POINT<TV,SCALAR>::TYPE> TYPE;};

template<class T> struct HAS_M_FOR_STATIC_SIZE{static const bool value=false;};
template<class T,int d> struct HAS_M_FOR_STATIC_SIZE<VECTOR<T,d> >{static const bool value=true;};
template<class T,int d> struct HAS_M_FOR_STATIC_SIZE<TWIST<VECTOR<T,d> > >{static const bool value=true;};

template<class T> struct IS_VECTOR{static const bool value=false;};
template<class T,int d> struct IS_VECTOR<ZERO_VECTOR<T,d> > {static const bool value=true;};
template<class T,int d> struct IS_VECTOR<VECTOR<T,d> > {static const bool value=true;};

}
#endif
