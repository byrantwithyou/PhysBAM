//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_FC_DRIVER__
#define __PLS_FC_DRIVER__
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class PLS_FC_EXAMPLE;

template<class TV>
class PLS_FC_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
public:
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    
    
    int current_frame;
    T time;
    int output_number;
    
    PLS_FC_EXAMPLE<TV>& example;
    
    PLS_FC_DRIVER(PLS_FC_EXAMPLE<TV>& example);
    virtual ~PLS_FC_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();
    
    void Advance_One_Time_Step(bool first_step);
    void Update_Pls(T dt);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Update_Level_Set(T dt,bool first_step);
    void Advection_And_BDF(T dt,bool first_step);
    void Apply_Pressure_And_Viscosity(T dt,bool first_step);
    void Extrapolate_Velocity(ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > >& u,const ARRAY<int,FACE_INDEX<TV::dimension> >& color);
    void Extrapolate_Velocity(ARRAY<T,FACE_INDEX<TV::dimension> >& u,const ARRAY<int,FACE_INDEX<TV::dimension> >& color,int c);
    void No_Advection_And_BDF(T dt,bool first_step,int c);
    void Reduced_Advection_And_BDF(T dt,bool first_step,int c);
    void RK2_Advection_And_BDF(T dt,bool first_step,int c);
    void Assert_Advection_CFL(const ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<int,FACE_INDEX<TV::m> >& color,int c,T dt) const;
//#####################################################################
};
}
#endif
