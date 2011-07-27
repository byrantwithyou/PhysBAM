//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONCONVEX_BOUNCE_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - January 4, 2003
//#####################################################################
#ifndef __NONCONVEX_BOUNCE_EXAMPLE__
#define __NONCONVEX_BOUNCE_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class NONCONVEX_BOUNCE_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    NONCONVEX_BOUNCE_EXAMPLE(int parameter=0);
    ~NONCONVEX_BOUNCE_EXAMPLE();
    void Initialize_Rigid_Bodies();
};
}
#endif
