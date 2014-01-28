//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_DRIVER
//#####################################################################
#ifndef __SOLIDS_DRIVER__
#define __SOLIDS_DRIVER__    

#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
namespace PhysBAM{

template<class TV>
class SOLIDS_DRIVER:public DRIVER<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef DRIVER<TV> BASE;
public:
    using BASE::output_number;using BASE::time;using BASE::Read_Time;using BASE::Write_Substep;
    using BASE::Write_Time;using BASE::Write_First_Frame;using BASE::Write_Last_Frame;
    SOLIDS_EXAMPLE<TV>& example;
    bool project_at_frame_boundaries;
    int current_frame;
    T next_dt; // for fluid time stepping
    bool next_done; // for fluid time stepping
    T last_dt;
    T restart_dt;
    bool reset_with_restart;

    SOLIDS_DRIVER(SOLIDS_EXAMPLE<TV>& example_input);
    virtual ~SOLIDS_DRIVER();

//#####################################################################
    virtual void Preprocess_Frame(const int frame);
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame_input) PHYSBAM_OVERRIDE;
    void Initialize() PHYSBAM_OVERRIDE;
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    T Compute_Dt(const T time,const T target_time,bool& done);
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Setup_Solids(const T time,const int substep);
    void Solid_Position_Update(const T dt,const int substep);
    void Rigid_Cluster_Fracture(const T dt_full_advance,const T dt_cfl,const int substep);
    void Solid_Velocity_Update(const T dt,const int substep,const bool done);
//#####################################################################
};
}
#endif
