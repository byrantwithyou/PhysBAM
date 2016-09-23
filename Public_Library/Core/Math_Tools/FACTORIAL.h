//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACTORIAL
//#####################################################################
#ifndef __FACTORIAL__
#define __FACTORIAL__

namespace PhysBAM{
constexpr inline int factorial(int d) {return d==0?1:d*factorial(d-1);}
}
#endif
