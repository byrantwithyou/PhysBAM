//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_FC_DRIVER__
#define __PLS_FC_DRIVER__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class PLS_FC_EXAMPLE;
template<class TV> class INTERFACE_STOKES_SYSTEM_COLOR;
template<class T> class KRYLOV_VECTOR_BASE;
template<class T,int d> class SYMMETRIC_MATRIX;

template<class TV>
class PLS_FC_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    
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
    void Update_Polymer_Stress(T dt);
    void Extrapolate_Velocity(ARRAY<ARRAY<T,FACE_INDEX<TV::m> > >& u,const ARRAY<int,FACE_INDEX<TV::m> >& color);
    void Extrapolate_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<int,FACE_INDEX<TV::m> >& color,int c);
    void Extrapolate_Polymer_Stress(ARRAY<ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> >& u,const ARRAY<int,TV_INT>& color);
    void Extrapolate_Polymer_Stress(ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& u,const ARRAY<int,TV_INT>& color,int c);
    void No_Advection_And_BDF(T dt,bool first_step,int c);
    void Reduced_Advection_And_BDF(T dt,bool first_step,int c);
    void RK2_Advection_And_BDF(T dt,bool first_step,int c);
    void Assert_Advection_CFL(const ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<int,FACE_INDEX<TV::m> >& color,int c,T dt) const;
    void Dump_Largest_Eigenvector(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,ARRAY<KRYLOV_VECTOR_BASE<T>*>& vectors) const;
    void Dump_Vector(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,const KRYLOV_VECTOR_BASE<T>& u,const char* str) const;
//#####################################################################
};
}
#endif
