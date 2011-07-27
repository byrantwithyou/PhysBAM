//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEST_EXAMPLE
//##################################################################### 
//
//#####################################################################
// Guendelman - May 2, 2003
//#####################################################################
#ifndef __TEST_EXAMPLE__
#define __TEST_EXAMPLE__

#include "../RIGID_BODIES_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class TEST_EXAMPLE:public RIGID_BODIES_EXAMPLE<T>
{
public:
    TEST_EXAMPLE(int parameter=0);
    ~TEST_EXAMPLE();
    void Initialize_Rigid_Bodies() PHYSBAM_OVERRIDE;

private:
    void Dribble_Test();
    void Chute_Test();
};
}
#endif
