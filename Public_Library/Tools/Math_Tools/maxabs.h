//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function maxabs
//#####################################################################
//
// finds the maximum absolute value
//
//#####################################################################
#ifndef __maxabs__
#define __maxabs__

#include <Tools/Math_Tools/max.h>
#include <cmath>
namespace PhysBAM{

using ::std::abs;

template<class T>
inline T maxabs(const T a,const T b)
{return max(abs(a),abs(b));}

template<class T,class ...Args>
inline T maxabs(const T a,const T b,Args&&... c)
{return max(abs(a),maxabs(b,c...));}

}
#endif
