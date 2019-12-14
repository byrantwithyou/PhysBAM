//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_DRIVER__
#define __PLS_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_PDE/Advection/ADVECTION_UNIFORM_FORWARD.h>
namespace PhysBAM{


template<class TV> class PLS_EXAMPLE;

template<class TV>
class PLS_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T> T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;

protected:
    int current_frame;
    T time;
    int output_number;

    PLS_EXAMPLE<TV>& example;
public:

    PLS_DRIVER(PLS_EXAMPLE<TV>& example);
    virtual ~PLS_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files();
    void Write_Substep(const std::string& title);

//#####################################################################
};
}
#endif
