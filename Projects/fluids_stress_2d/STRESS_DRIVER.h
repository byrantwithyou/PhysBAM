//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STRESS_DRIVER__
#define __STRESS_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
namespace PhysBAM{

template<class TV> class STRESS_EXAMPLE;
template<class TV> class INTERFACE_STOKES_SYSTEM_COLOR;
template<class T> class KRYLOV_VECTOR_BASE;

template<class TV>
class STRESS_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    
    int current_frame;
    T time;
    int output_number;
    
    STRESS_EXAMPLE<TV>& example;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> next_polymer_stress;
    ARRAY<T,FACE_INDEX<TV::m> > temp_face_velocities,temp_face_velocities2;

    STRESS_DRIVER(STRESS_EXAMPLE<TV>& example);
    virtual ~STRESS_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();
    
    void Advance_One_Time_Step(bool first_step);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Advection(T dt,bool one_step,int from_time,int to_time,int bc_time);
    void Extrapolate_Stress(ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& S);
    void Assert_Advection_CFL(const ARRAY<T,FACE_INDEX<TV::m> >& u,T dt) const;
    void Add_RHS_Terms(T dt);
//#####################################################################
};
}
#endif
