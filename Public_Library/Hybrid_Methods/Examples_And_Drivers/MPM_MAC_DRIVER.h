//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MAC_DRIVER__
#define __MPM_MAC_DRIVER__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <climits>
#include <functional>
namespace PhysBAM{

template<class TV> class MPM_MAC_EXAMPLE;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class LAPLACE_UNIFORM;
enum class RF;

template<class TV>
class MPM_MAC_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
public:
    enum {pressure_D=-1,pressure_N=-2,pressure_uninit=-3};

    int current_frame;
    int output_number;

    MPM_MAC_EXAMPLE<TV>& example;

    MPM_MAC_DRIVER(MPM_MAC_EXAMPLE<TV>& example);
    virtual ~MPM_MAC_DRIVER();

    void Execute_Main_Program();
    void Initialize();

    void Advance_One_Time_Step();
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files();
    void Write_Substep(const std::string& title);
    void Update_Particle_Weights();
    void Prepare_Scatter();
    void Particle_To_Grid();
    void Grid_To_Particle();
    void Compute_Effective_Velocity();
    void Build_Level_Sets();
    void Pressure_Projection();
    void Apply_Forces();
    void Apply_Viscosity();
    void Move_Particles();
    void Extrapolate_Velocity(bool use_bc,bool extrapolate_boundary);
    template <class D,class N>
    void Reflect_Boundary(D func_d,N func_n,RF flag) const;
    void Reflect_Boundary_Mass_Momentum(ARRAY<T,FACE_INDEX<TV::m> >& P) const;
    void Reflect_Boundary_Particle_Force(ARRAY<T,FACE_INDEX<TV::m> >& force) const;
    void Reflect_Boundary_Grid_Force(ARRAY<T,FACE_INDEX<TV::m> >& force) const;
    void Reflect_Boundary_Velocity_Copy_Only() const;
    void Extrapolate_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& velocity,bool use_bc,bool extrapolate_boundary) const;
    void Extrapolate_Boundary(ARRAY<T,FACE_INDEX<TV::m> >& velocity) const;
    T Compute_Dt() const;
    T Max_Particle_Speed() const;
    T Grid_V_Upper_Bound() const;
    void Update_Simulated_Particles();
    void Compute_Poisson_Matrix();
    template<class T2> void Fix_Periodic(ARRAY<T2,TV_INT>& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic_Node(ARRAY<T2,TV_INT>& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic_Accum(ARRAY<T2,TV_INT>& u,int ghost=INT_MAX) const;
    template<class T2> void Fix_Periodic_Accum(ARRAY<T2,FACE_INDEX<TV::m> >& u,int ghost=INT_MAX) const;
    T Density(const FACE_INDEX<TV::m>& face_index) const;

    void Compute_Boundary_Conditions();
    int Allocate_Projection_System_Variable();
    void Compute_Laplacian(int var);
    void Compute_Gradient(int nvar);
    void Reseeding();
    void Step(std::function<void()> func,const char* name,bool dump_substep=true,bool do_step=true);
    void Invalidate_Particle(int p);
    bool Neumann_Boundary_Condition(const FACE_INDEX<TV::m>& face,T& bc) const;
    void Extrapolate_Inside_Object();
    void Level_Set_Pressure_Projection();
//#####################################################################
};
}
#endif
