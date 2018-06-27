//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STOKES_MF_DRIVER__
#define __STOKES_MF_DRIVER__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class STOKES_MF_EXAMPLE;

template<class TV>
class STOKES_MF_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    enum {pressure_D=-1,pressure_N=-2,pressure_uninit=-3};

    int current_frame;
    int output_number;

    STOKES_MF_EXAMPLE<TV>& example;

    STOKES_MF_DRIVER(STOKES_MF_EXAMPLE<TV>& example);
    virtual ~STOKES_MF_DRIVER();

    void Execute_Main_Program();
    void Initialize();
//#####################################################################
};
}
#endif
