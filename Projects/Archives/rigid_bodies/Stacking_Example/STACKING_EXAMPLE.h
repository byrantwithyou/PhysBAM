//#####################################################################
// Copyright 2002,2003, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STACKING_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - January 1, 2003
//#####################################################################
#ifndef __STACKING_EXAMPLE__
#define __STACKING_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class STACKING_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    STACKING_EXAMPLE(int parameter=0);
    ~STACKING_EXAMPLE();
    void Initialize_Rigid_Bodies();
};
}
#endif
