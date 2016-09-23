//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ARRAYS_UNIFORM_FORWARD
//#####################################################################
#ifndef __ARRAYS_UNIFORM_FORWARD__
#define __ARRAYS_UNIFORM_FORWARD__

#include <Core/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> class ARRAY_BASE;
template<class T,class TV_INT> class ARRAYS_ND_BASE;
template<class T,class ID> class ARRAY;
template<int d> class FACE_INDEX;
template<int d> class SIDED_FACE_INDEX;

template<int d> class FLOOD_FILL;

template<class T,int d> class VECTOR;

}
#endif
