//#####################################################################
// Copyright 2002, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SLIDING_TEST
//##################################################################### 
//
//#####################################################################
// Guendelman - December 3, 2002
//#####################################################################
#ifndef __SLIDING_TEST__
#define __SLIDING_TEST__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T> 
class SLIDING_TEST:public RIGID_BODIES_EXAMPLE<T>
{
public:
    SLIDING_TEST(int parameter=0);
    ~SLIDING_TEST();
    void Initialize_Rigid_Bodies();
};
}
#endif
