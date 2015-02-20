//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_EXAMPLE__
#define __MPM_EXAMPLE__
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class PARTICLE_GRID_FORCES;
template<class T> class KRYLOV_VECTOR_BASE;
template<class T> class MPM_KRYLOV_VECTOR;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class GATHER_SCATTER;
template<class TV> class DEBUG_PARTICLES;
template<class TV> class IMPLICIT_OBJECT;
template<class TV> class MPM_COLLISION_OBJECT;
template<class TV> class DEFORMABLES_FORCES;

template<class TV>
class MPM_EXAMPLE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    GRID<TV> grid;
    STREAM_TYPE stream_type;
    MPM_PARTICLES<TV>& particles;
    DEBUG_PARTICLES<TV>& debug_particles;
    ARRAY<int> simulated_particles;

    ARRAY<T,TV_INT> mass;
    ARRAY<TV,TV_INT> velocity,velocity_new;
    ARRAY<int> valid_grid_indices;
    ARRAY<PARTICLE_GRID_FORCES<TV>*> forces;
    ARRAY<DEFORMABLES_FORCES<TV>*> lagrangian_forces;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    MPM_KRYLOV_VECTOR<TV>& rhs;
    PARTICLE_GRID_WEIGHTS<TV>* weights;
    GATHER_SCATTER<TV>& gather_scatter;
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    mutable ARRAY<TV> lagrangian_forces_V,lagrangian_forces_F;

    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    std::string output_directory;
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    int ghost;
    bool use_reduced_rasterization;
    bool use_affine;
    bool use_midpoint;
    bool use_particle_collision;
    T flip;
    T cfl;

    T newton_tolerance;
    int newton_iterations;
    T solver_tolerance;
    int solver_iterations;
    bool test_diff;
    int threads;

    MPM_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~MPM_EXAMPLE();
    
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    virtual void Begin_Frame(const int frame)=0;
    virtual void End_Frame(const int frame)=0;
    virtual void Begin_Time_Step(const T time)=0;
    virtual void End_Time_Step(const T time)=0;

    void Precompute_Forces(const T time);
    T Potential_Energy(const T time) const;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const;
    void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const;
    int Add_Force(PARTICLE_GRID_FORCES<TV>& force);
    int Add_Force(DEFORMABLES_FORCES<TV>& force);
    void Set_Weights(PARTICLE_GRID_WEIGHTS<TV>* weights_input);

//#####################################################################
};
}
#endif
