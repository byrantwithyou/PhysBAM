//#####################################################################
// Copyright 2002, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNCE_FRICTION_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - December 18, 2002
//#####################################################################
#ifndef __BOUNCE_FRICTION_EXAMPLE__
#define __BOUNCE_FRICTION_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class BOUNCE_FRICTION_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    BOUNCE_FRICTION_EXAMPLE(int parameter=0);
    ~BOUNCE_FRICTION_EXAMPLE();
    void Initialize_Rigid_Bodies();
};
}
#endif
