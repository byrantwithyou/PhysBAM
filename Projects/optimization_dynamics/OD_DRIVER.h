//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OD_DRIVER__
#define __OD_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include "OD_EXAMPLE.h"
#include "OD_SOLVER.h"
namespace PhysBAM{

template<class TV>
class OD_DRIVER
{
    typedef typename TV::SCALAR T;
public:

    int current_frame;
    int output_number;

    OD_EXAMPLE<TV>& example;
    OD_SOLVER<TV>& solver;

    OD_DRIVER(OD_EXAMPLE<TV>& example,OD_SOLVER<TV>& solver);
    virtual ~OD_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Apply_External_Forces();
//#####################################################################
};
}

#endif
