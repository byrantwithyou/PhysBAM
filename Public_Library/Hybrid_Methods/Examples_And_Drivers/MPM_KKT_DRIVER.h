//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_KKT_DRIVER__
#define __MPM_KKT_DRIVER__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV> class MPM_KKT_KRYLOV_SYSTEM;
template<class TV> class MPM_KKT_KRYLOV_VECTOR;
template<class TV> class MPM_KKT_EXAMPLE;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class T> class KRYLOV_VECTOR_BASE;
template<class T> class GRID;

template<class TV>
class MPM_KKT_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
public:

    int current_frame;
    int output_number;

    MPM_KKT_EXAMPLE<TV>& example;
    MPM_KKT_KRYLOV_SYSTEM<TV>& kkt_sys;
    MPM_KKT_KRYLOV_VECTOR<TV>& kkt_lhs;
    MPM_KKT_KRYLOV_VECTOR<TV>& kkt_rhs;
    
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;

    MPM_KKT_DRIVER(MPM_KKT_EXAMPLE<TV>& example);
    virtual ~MPM_KKT_DRIVER();
    
    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Update_Particle_Weights();
    void Particle_To_Grid();
    void Grid_To_Particle();
    void Solve_KKT_System();
    void Apply_Forces();
    void Apply_Friction();
    void Initialize_Location(const RANGE<TV_INT>& range,const GRID<TV>& grid,ARRAY<TV,TV_INT>& location);
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Grid_V_Upper_Bound() const;
    void Update_Simulated_Particles();
    void Print_Grid_Stats(const char* str,T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0);
    void Print_Particle_Stats(const char* str,T dt);
    void Print_Energy_Stats(const char* str,const ARRAY<TV,TV_INT>& u);
    void Build_FEM_Mass_Matrix();
//#####################################################################
};
}
#endif
