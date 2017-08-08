//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_EXAMPLE__
#define __MPM_EXAMPLE__
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <functional>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class GATHER_SCATTER;
template<class TV> class IMPLICIT_OBJECT;
template<class T>  class KRYLOV_VECTOR_BASE;
template<class TV> class MPM_FORCE_HELPER;
template<class T>  class MPM_KRYLOV_VECTOR;
template<class TV> class MPM_PARTICLES;
template<class TV> class PARTICLE_GRID_FORCES;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_PLASTICITY_MODEL;

template<class TV>
class MPM_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
    GRID<TV> grid;
    STREAM_TYPE stream_type;
    MPM_PARTICLES<TV>& particles;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    DEBUG_PARTICLES<TV>& debug_particles;
    ARRAY<int> simulated_particles;
    ARRAY<bool> particle_is_simulated;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > Dp_inv;
    ARRAY<T,TV_INT> mass,volume;
    ARRAY<TV,TV_INT> location;
    ARRAY<TV,TV_INT> velocity,velocity_new,velocity_friction,*current_velocity;
    ARRAY<int> valid_grid_indices;
    ARRAY<TV_INT> valid_grid_cell_indices;
    ARRAY<PARTICLE_GRID_FORCES<TV>*> forces;
    ARRAY<MPM_PLASTICITY_MODEL<TV>*> plasticity_models;
    ARRAY<DEFORMABLES_FORCES<TV>*>& lagrangian_forces;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    PARTICLE_GRID_WEIGHTS<TV>* weights;
    GATHER_SCATTER<TV>& gather_scatter;
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    mutable ARRAY<TV> lagrangian_forces_V,lagrangian_forces_F;
    MPM_FORCE_HELPER<TV>& force_helper;

    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    std::string output_directory;
    std::string data_directory;
    std::string test_output_prefix;
    bool use_test_output;
    T mass_contour;
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    int ghost;
    bool use_affine;
    bool use_midpoint;
    bool use_symplectic_euler;
    bool lag_Dp; // Bp will actually store Cp
    bool use_oldroyd;
    bool print_stats;
    bool only_write_particles;
    T flip;
    T cfl;
    T inv_Wi;

    T newton_tolerance;
    int newton_iterations;
    T solver_tolerance;
    int solver_iterations;
    bool test_diff;
    int threads;

    TV last_linear_momentum;
    typename TV::SPIN last_angular_momentum;
    T last_te;
    T last_grid_ke;
    bool output_structures_each_frame;
    T quad_F_coeff;
    bool asymmetric_system;

    MPM_EXAMPLE(const STREAM_TYPE stream_type_input);
    MPM_EXAMPLE(const MPM_EXAMPLE&) = delete;
    void operator=(const MPM_EXAMPLE&) = delete;
    virtual ~MPM_EXAMPLE();
    
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    std::function<void(int frame)> begin_frame;
    std::function<void(int frame)> end_frame;
    std::function<void(T time)> begin_time_step;
    std::function<void(T time)> end_time_step;

    void Capture_Stress();
    void Precompute_Forces(const T time,const T dt,const bool update_hessian);
    T Potential_Energy(const T time) const;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const;
    void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const;
    void Update_Lagged_Forces(const T time) const;
    int Add_Force(PARTICLE_GRID_FORCES<TV>& force);
    int Add_Force(DEFORMABLES_FORCES<TV>& force);
    void Set_Weights(int order);
    void Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0);
    template<class OBJECT> typename enable_if<!is_pointer<OBJECT>::value>::type
    Add_Collision_Object(const OBJECT& object,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0)
    {Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<OBJECT>(object),type,friction,func_frame,func_twist);}

    TV Total_Particle_Linear_Momentum() const;
    TV Total_Grid_Linear_Momentum(const ARRAY<TV,TV_INT>& u) const;
    typename TV::SPIN Total_Grid_Angular_Momentum(T dt,const ARRAY<TV,TV_INT>& u,const ARRAY<TV,TV_INT>* u0) const;
    typename TV::SPIN Total_Particle_Angular_Momentum() const;
    T Total_Grid_Kinetic_Energy(const ARRAY<TV,TV_INT>& u) const;
    T Total_Particle_Kinetic_Energy() const;
    T Average_Particle_Mass() const;
//#####################################################################
};
}
#endif
