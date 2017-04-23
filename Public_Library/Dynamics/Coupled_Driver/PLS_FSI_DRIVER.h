//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLS_FSI_DRIVER
//#####################################################################
#ifndef __PLS_FSI_DRIVER__
#define __PLS_FSI_DRIVER__    

#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
namespace PhysBAM{

template<class TV>
class PLS_FSI_DRIVER:public DRIVER<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef DRIVER<TV> BASE;
    typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    using BASE::time;
public:
    using BASE::output_number;using BASE::Write_Output_Files;using BASE::Read_Time;using BASE::Write_Time;
    using BASE::Write_First_Frame;using BASE::Write_Last_Frame;using BASE::Write_Substep;

    PLS_FSI_EXAMPLE<TV>& example;
    int current_frame;
    ARRAY<T,TV_INT> old_phi;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_scratch;

    PLS_FSI_DRIVER(PLS_FSI_EXAMPLE<TV>& example_input);
    PLS_FSI_DRIVER(const PLS_FSI_DRIVER&) = delete;
    void operator=(const PLS_FSI_DRIVER&) = delete;
    virtual ~PLS_FSI_DRIVER();

//#####################################################################
    void Initialize() override;
    void Initialize_Fluids_Grids();
    void Advance_To_Target_Time(const T target_time) override;
    void First_Order_Time_Step(int substep,T dt);
    virtual void Postprocess_Frame(const int frame);
    virtual void Preprocess_Frame(const int frame);
    T Compute_Dt(const T time,const T target_time,bool& done);
    void Write_Output_Files(const int frame) override;
    void Advect_Fluid(const T dt,const int substep);
    void Execute_Main_Program() override;
    void Simulate_To_Frame(const int frame_input) override;
    void Delete_Particles_Inside_Objects(const T time);
    template<class T_PARTICLES> void Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    void Extrapolate_Velocity_Across_Interface(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const LEVELSET<TV>& phi,const T band_width);
    void Advance_Particles_With_PLS(T dt);
    void Extrapolate_Velocity_Across_Interface(T time,T dt);
//#####################################################################
};
}
#endif
