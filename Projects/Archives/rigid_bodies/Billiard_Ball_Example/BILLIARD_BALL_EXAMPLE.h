//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BILLIARD_BALL_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - January 4, 2003
//#####################################################################
#ifndef __BILLIARD_BALL_EXAMPLE__
#define __BILLIARD_BALL_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class BILLIARD_BALL_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    BILLIARD_BALL_EXAMPLE(int parameter=0);
    ~BILLIARD_BALL_EXAMPLE();
    void Initialize_Rigid_Bodies();
};
}
#endif
