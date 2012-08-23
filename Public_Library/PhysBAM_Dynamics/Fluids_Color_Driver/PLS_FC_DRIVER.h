//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_FC_DRIVER__
#define __PLS_FC_DRIVER__
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
namespace PhysBAM{


template<class TV> class PLS_FC_EXAMPLE;

template<class TV>
class PLS_FC_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;

protected:
    int current_frame;
    T time;
    int output_number;

    PLS_FC_EXAMPLE<TV>& example;
public:

    PLS_FC_DRIVER(PLS_FC_EXAMPLE<TV>& example);
    virtual ~PLS_FC_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step(bool first_step);
    void Update_Pls(T dt);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Advection_And_BDF(T dt,bool first_step);
    void Apply_Pressure_And_Viscosity(T dt,bool first_step);

//#####################################################################
};
}
#endif
