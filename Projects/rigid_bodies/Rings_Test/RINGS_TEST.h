//#####################################################################
// Copyright 2002, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RINGS_TEST
//##################################################################### 
//
//#####################################################################
// Guendelman - December 3, 2002
//#####################################################################
#ifndef __RINGS_TEST__
#define __RINGS_TEST__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class RINGS_TEST:public RIGID_BODIES_EXAMPLE<T>
{
public:
    RINGS_TEST(int parameter=0);
    ~RINGS_TEST();
    void Initialize_Rigid_Bodies();
};
}
#endif
