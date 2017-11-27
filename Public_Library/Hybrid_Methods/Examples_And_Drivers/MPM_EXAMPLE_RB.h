//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_EXAMPLE_RB__
#define __MPM_EXAMPLE_RB__
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
namespace PhysBAM{

template<class TV>
class MPM_EXAMPLE_RB:public MPM_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_EXAMPLE_RB(const STREAM_TYPE stream_type_input);
    MPM_EXAMPLE_RB(const MPM_EXAMPLE_RB&) = delete;
    void operator=(const MPM_EXAMPLE_RB&) = delete;
    virtual ~MPM_EXAMPLE_RB();
//#####################################################################
};
}
#endif
