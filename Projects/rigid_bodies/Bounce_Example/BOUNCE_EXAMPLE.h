//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNCE_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - January 2, 2003
//#####################################################################
#ifndef __BOUNCE_EXAMPLE__
#define __BOUNCE_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class BOUNCE_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    BOUNCE_EXAMPLE(int parameter=0);
    ~BOUNCE_EXAMPLE();
    void Initialize_Rigid_Bodies();
};
}
#endif
