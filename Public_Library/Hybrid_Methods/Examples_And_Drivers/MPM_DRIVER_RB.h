//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_DRIVER_RB__
#define __MPM_DRIVER_RB__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV> class MPM_EXAMPLE_RB;
template<class TV> class MPM_OBJECTIVE_RB;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_KRYLOV_VECTOR_RB;
template<class T> class KRYLOV_VECTOR_BASE;

template<class TV>
class MPM_DRIVER_RB
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
public:

    int current_frame;
    int output_number;

    MPM_EXAMPLE_RB<TV>& example;
    MPM_OBJECTIVE_RB<TV>& objective;
    MPM_KRYLOV_VECTOR_RB<TV>& dv,&rhs;
    ARRAY<TWIST<TV> > rigid_forces;
    
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    
    MPM_DRIVER_RB(MPM_EXAMPLE_RB<TV>& example);
    virtual ~MPM_DRIVER_RB();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title);
    void Update_Particle_Weights();
    void Particle_To_Grid();
    void Register_Particles();
    void Grid_To_Particle();
    void Update_Plasticity_And_Hardening();
    void Apply_Forces();
    void Apply_Friction();
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Grid_V_Upper_Bound() const;
    void Update_Simulated_Particles();
    void Print_Grid_Stats(const char* str,T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0);
    void Print_Particle_Stats(const char* str,T dt);
    void Print_Energy_Stats(const char* str,const ARRAY<TV,TV_INT>& u);
    void Grid_To_Particle_Limit_Dt();
    void Limit_Dt_Sound_Speed();
    void Reflect_Or_Invalidate_Particle(int p);
    void Move_Rigid_Bodies();
    void Apply_Rigid_Body_Forces();
    void Process_Pairwise_Collisions();
    void Process_Projected_Collisions(T dt);
//#####################################################################
};

template<class TV,class T>
FRAME<TV> Move_Rigid_Body(T dt,const FRAME<TV>& frame,const TWIST<TV>& twist,
    const SYMMETRIC_MATRIX<T,TV::SPIN::m>& inertia);
}
#endif
