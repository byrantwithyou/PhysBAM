//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#ifndef __COMMON_DATA__
#define __COMMON_DATA__
#include <Core/Vectors/VECTOR.h>

namespace PhysBAM
{
template<class TV>
struct COMMON_DATA
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    
    int test_number=-1,resolution=32,threads=1;
    T m=1,s=1,kg=1;
    bool user_resolution=false;
    bool user_last_frame=false;
    int seed=1234;
    
    COMMON_DATA()
    {
    }
};
}

#endif
