//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MICROPOLAR_DRIVER__
#define __MPM_MICROPOLAR_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV> class FLUID_KRYLOV_SYSTEM;
template<class TV> class FLUID_KRYLOV_VECTOR;
template<class TV> class MPM_MICROPOLAR_EXAMPLE;
template<class TV> class MPM_OBJECTIVE;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_KRYLOV_VECTOR;
template<class T> class KRYLOV_VECTOR_BASE;

template<class TV>
class MPM_MICROPOLAR_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
public:

    int current_frame;
    int output_number;

    MPM_MICROPOLAR_EXAMPLE<TV>& example;
    ARRAY<TV,TV_INT> forces;
    
    MPM_MICROPOLAR_DRIVER(MPM_MICROPOLAR_EXAMPLE<TV>& example);
    virtual ~MPM_MICROPOLAR_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title);
    void Update_Particle_Weights();
    void Particle_To_Grid();
    void Grid_To_Particle();
    void Apply_Forces();
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Grid_V_Upper_Bound() const;
    void Update_Simulated_Particles();
    void Grid_To_Particle_Limit_Dt();
    void Limit_Dt_Sound_Speed();
    template<class S> void Reflection_Boundary_Condition(ARRAY<S,TV_INT>& u,bool flip_sign);
    void Reflect_Or_Invalidate_Particle(int p);
//#####################################################################
};
}
#endif
