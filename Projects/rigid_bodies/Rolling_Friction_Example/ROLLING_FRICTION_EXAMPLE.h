//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROLLING_FRICTION_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - January 5, 2003
//#####################################################################
#ifndef __ROLLING_FRICTION_EXAMPLE__
#define __ROLLING_FRICTION_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class ROLLING_FRICTION_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    ROLLING_FRICTION_EXAMPLE(int parameter=0);
    ~ROLLING_FRICTION_EXAMPLE();
    void Initialize_Rigid_Bodies();
};
}
#endif
