//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLANK_TEST
//##################################################################### 
//
//#####################################################################
// Guendelman - January 2, 2003
//#####################################################################
#ifndef __PLANK_TEST__
#define __PLANK_TEST__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class PLANK_TEST:public RIGID_BODIES_EXAMPLE<T>
{
public:
    PLANK_TEST(int parameter=0);
    ~PLANK_TEST();
    void Initialize_Rigid_Bodies();
};
}
#endif
