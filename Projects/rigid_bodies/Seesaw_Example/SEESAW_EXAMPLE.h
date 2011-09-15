//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEESAW_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - May 2, 2003
//#####################################################################
#ifndef __SEESAW_EXAMPLE__
#define __SEESAW_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T> 
class SEESAW_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    SEESAW_EXAMPLE(int parameter=0);
    ~SEESAW_EXAMPLE();
    void Initialize_Rigid_Bodies();
};
}
#endif
