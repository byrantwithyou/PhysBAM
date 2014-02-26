//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function minabs
//#####################################################################
//
// finds the minimum absolute value
//
//#####################################################################
#ifndef __minabs__
#define __minabs__

#include <Tools/Math_Tools/min.h>
#include <cmath>
namespace PhysBAM{

template<class T>
inline T minabs(const T a,const T b)
{return min(abs(a),abs(b));}

template<class T,class ...Args>
inline T minabs(const T a,const T b,Args&&... c)
{return min(abs(a),minabs(b,c...));}

}
#endif
