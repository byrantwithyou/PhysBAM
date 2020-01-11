//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_DRIVER__
#define __MPM_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV> class FLUID_KRYLOV_SYSTEM;
template<class TV> class FLUID_KRYLOV_VECTOR;
template<class TV> class MPM_EXAMPLE;
template<class TV> class MPM_OBJECTIVE;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_KRYLOV_VECTOR;
template<class T> class KRYLOV_VECTOR_BASE;
template<class TV> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;
enum class RI;

template<class TV>
class MPM_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
public:

    int current_frame=0;

    MPM_EXAMPLE<TV>& example;
    MPM_OBJECTIVE<TV>& objective;
    MPM_KRYLOV_VECTOR<TV>& dv,&rhs;
    
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    MPM_DRIVER(MPM_EXAMPLE<TV>& example);
    virtual ~MPM_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files();
    void Write_Substep(const std::string& title);
    void Update_Particle_Weights();
    void Particle_To_Grid();
    void Grid_To_Particle();
    void Update_Plasticity_And_Hardening();
    void Apply_Forces();
    void Apply_Particle_Forces(ARRAY<TV,TV_INT>& F);
    void Apply_Grid_Forces(ARRAY<TV,TV_INT>& F);
    void Apply_Friction();
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Grid_V_Upper_Bound() const;
    T Compute_Max_Sound_Speed() const;
    T Max_Dt_Single_Particle_Pressure() const;
    T Max_Dt_Single_Particle() const;
    void Update_Simulated_Particles();
    void Print_Grid_Stats(const char* str,T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0);
    void Print_Particle_Stats(const char* str,T dt);
    void Print_Energy_Stats(const char* str,const ARRAY<TV,TV_INT>& u);
    T Grid_To_Particle_Limit_Dt();
    T Limit_Dt_Sound_Speed();
    void Reduce_Dt();
    void Reflect_Boundary_Mass_Momentum();
    void Reflect_Boundary_Friction(const ARRAY<T,TV_INT>& mass,ARRAY<TV,TV_INT>& u) const;
    void Reflect_Boundary_Velocity(ARRAY<TV,TV_INT>& u);
    void Reflect_Boundary_Force(ARRAY<TV,TV_INT>& force);
    template<class S, class L, class N> void Reflect_Boundary(S func_s,L func_l,N func_n,RI flag) const;
    void Reflect_Or_Invalidate_Particle(int p);
    void Test_Sound_Speed(int num) const;
    void Step(std::function<void()> func,const char* name,bool dump_substep=true,bool do_step=true);
    void Apply_Reflection_Collision_Objects();
    void Sample_Reflection_Collision_Object(int i);
//#####################################################################
};
}
#endif
