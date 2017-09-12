//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SHALLOW_WATER_DRIVER__
#define __SHALLOW_WATER_DRIVER__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class SHALLOW_WATER_STATE;

template<class TV>
class SHALLOW_WATER_DRIVER
{
    STATIC_ASSERT(TV::m<3);
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:

    int current_frame;
    int output_number;

    SHALLOW_WATER_STATE<TV>& state;
    
    SHALLOW_WATER_DRIVER(SHALLOW_WATER_STATE<TV>& state);
    virtual ~SHALLOW_WATER_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title);
    void Apply_Forces();
    T Compute_Dt() const;
//#####################################################################
};
}
#endif
