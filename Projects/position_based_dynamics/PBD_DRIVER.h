//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PBD_DRIVER__
#define __PBD_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include "PBD_EXAMPLE.h"
namespace PhysBAM{

template<class TV>
class PBD_DRIVER
{
    typedef typename TV::SCALAR T;
public:

    int current_frame;
    int output_number;

    PBD_EXAMPLE<TV>& example;

    PBD_DRIVER(PBD_EXAMPLE<TV>& example);
    virtual ~PBD_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title);
    void Apply_External_Forces();
    void Project_Constraints();
//#####################################################################
};
}

#endif
