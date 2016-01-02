//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_KKT_EXAMPLE__
#define __MPM_KKT_EXAMPLE__
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;
template<class TV> class PARTICLE_GRID_FORCES;
template<class T> class KRYLOV_VECTOR_BASE;
template<class T> class KKT_KRYLOV_VECTOR;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class GATHER_SCATTER;
template<class TV> class DEBUG_PARTICLES;
template<class TV> class IMPLICIT_OBJECT;
template<class TV> class MPM_COLLISION_OBJECT;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class MPM_COLLISION_OBJECT;
template<class TV> class MPM_FORCE_HELPER;

template<class TV>
class MPM_KKT_EXAMPLE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
    GRID<TV> grid;
    GRID<TV> coarse_grid;
    STREAM_TYPE stream_type;
    MPM_PARTICLES<TV>& particles;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    DEBUG_PARTICLES<TV>& debug_particles;
    ARRAY<int> simulated_particles;
    ARRAY<bool> particle_is_simulated;
    ARRAY<T,TV_INT> mass,mass_coarse;
    ARRAY<TV,TV_INT> location;
    ARRAY<TV,TV_INT> coarse_location;
    ARRAY<TV,TV_INT> velocity,velocity_new,*current_velocity;
    ARRAY<int> valid_velocity_indices;
    ARRAY<TV_INT> valid_velocity_cell_indices;
    ARRAY<int> valid_pressure_indices;
    ARRAY<TV_INT> valid_pressure_cell_indices;
    ARRAY<int,TV_INT> inv_valid_pressure_cell;
    ARRAY<int,TV_INT> inv_valid_velocity_cell;
    ARRAY<PARTICLE_GRID_FORCES<TV>*> forces;
    ARRAY<DEFORMABLES_FORCES<TV>*>& lagrangian_forces;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    PARTICLE_GRID_WEIGHTS<TV>* weights;
    PARTICLE_GRID_WEIGHTS<TV>* coarse_weights;
    GATHER_SCATTER<TV>& gather_scatter;
    GATHER_SCATTER<TV>& gather_scatter_coarse;
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    ARRAY<IMPLICIT_OBJECT<TV>* > fluid_walls;
    mutable ARRAY<TV> lagrangian_forces_V,lagrangian_forces_F;
    MPM_FORCE_HELPER<TV>& force_helper;

    // kkt stuff
    SPARSE_MATRIX_FLAT_MXN<T> M;
    ARRAY<int,TV_INT> cell_flags;
    ARRAY<T,TV_INT> one_over_lambda;
    ARRAY<T,TV_INT> J;

    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    std::string output_directory;
    std::string data_directory;
    T mass_contour;
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    int ghost;
    bool use_reduced_rasterization;
    bool use_affine;
    bool use_f2p;
    bool use_midpoint;
    bool use_symplectic_euler;
    bool use_early_gradient_transfer;
    bool use_oldroyd;
    bool print_stats;
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
    bool use_FEM_mass;

    MPM_KKT_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~MPM_KKT_EXAMPLE();
    
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    virtual void Begin_Frame(const int frame)=0;
    virtual void End_Frame(const int frame)=0;
    virtual void Begin_Time_Step(const T time)=0;
    virtual void End_Time_Step(const T time)=0;

    void Capture_Stress();
    void Precompute_Forces(const T time,const T dt,const bool update_hessian);
    T Potential_Energy(const T time) const;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const;
    void Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const;
    void Update_Lagged_Forces(const T time) const;
    int Add_Force(PARTICLE_GRID_FORCES<TV>& force);
    int Add_Force(DEFORMABLES_FORCES<TV>& force);
    void Set_Weights(PARTICLE_GRID_WEIGHTS<TV>* weights_input);
    void Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction);
    void Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io);
    template<class OBJECT> typename enable_if<!is_pointer<OBJECT>::value>::type
    Add_Collision_Object(const OBJECT& object,COLLISION_TYPE type,T friction)
    {Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<OBJECT>(object),type,friction);}

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
