//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_EXAMPLE__
#define __MPM_EXAMPLE__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <functional>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class SOLID_BODY_COLLECTION;
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
    enum BC_TYPE {BC_FREE, BC_SLIP, BC_NOSLIP, BC_SEP};
    GRID<TV> grid;
    STREAM_TYPE stream_type;
    MPM_PARTICLES<TV>& particles;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    DEBUG_PARTICLES<TV>& debug_particles;
    ARRAY<int> simulated_particles;
    ARRAY<bool> particle_is_simulated;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > Dp_inv;
    ARRAY<T,TV_INT> mass,volume;
    ARRAY<TV,TV_INT> location;
    ARRAY<TV,TV_INT> velocity,velocity_save,velocity_friction_save;
    ARRAY<int> valid_grid_indices;
    ARRAY<TV_INT> valid_grid_cell_indices;
    ARRAY<PARTICLE_GRID_FORCES<TV>*> forces;
    ARRAY<MPM_PLASTICITY_MODEL<TV>*> plasticity_models;
    ARRAY<DEFORMABLES_FORCES<TV>*>& lagrangian_forces;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    PARTICLE_GRID_WEIGHTS<TV>* weights=0;
    GATHER_SCATTER<TV>& gather_scatter;
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    RANDOM_NUMBERS<T> random;
    ARRAY<int,TV_INT> colliding_nodes;

    struct REFLECT_OBJECT_DATA
    {
        ARRAY<TV> X;
    };

    ARRAY<REFLECT_OBJECT_DATA*> collision_objects_reflection;
    mutable ARRAY<TV> lagrangian_forces_V,lagrangian_forces_F;
    MPM_FORCE_HELPER<TV>& force_helper;

    HASHTABLE<std::string,PAIR<bool,VECTOR<ARRAY<std::function<void()> >,2> > > time_step_callbacks; // begin, end

    int last_frame=100;
    std::string frame_title;
    int write_substeps_level=-1;
    int substeps_delay_frame=-1;
    VIEWER_DIR viewer_dir{"output"};
    std::string data_directory="../../Public_Data";
    std::string test_output_prefix;
    bool use_test_output=false;
    T mass_contour=-1;
    bool auto_restart=false;
    int restart=0;
    T dt=0,time=0,frame_dt=(T)1/24;
    T min_dt=0,max_dt=frame_dt;
    int ghost=3;
    bool use_affine=true;
    bool use_midpoint=false;
    bool use_symplectic_euler=false;
    bool lag_Dp=false; // Bp will actually store Cp
    bool use_oldroyd=false;
    bool print_stats=false;
    bool only_write_particles=false;
    T flip=0;
    T cfl=1;
    T cfl_F=(T).5;
    T cfl_sound=(T).9;
    T inv_Wi=0;
    bool r_sound_speed=false;
    bool r_cfl=false;
    bool r_F=false;
    bool extra_render=false;
    bool use_strong_cfl=false;
    bool use_sound_speed_cfl=false;
    bool compute_sound_speed=false;
    bool reflection_bc=false;
    int reflection_bc_flags=0;
    bool use_full_reflection=false;
    T reflection_bc_friction=0;
    bool use_reflection_collision_objects=false;
    bool test_sound_speed=false;
    bool dilation_only=false;
    bool write_structures_every_frame=false;
    bool use_single_particle_cfl=false;
    T cfl_single_particle=(T)0.9;
    bool use_affine_cfl=true;
    bool verbose_cfl=true;
    
    T newton_tolerance=1;
    int newton_iterations=100;
    T solver_tolerance=(T).5;
    int solver_iterations=1000;
    bool test_diff=false;
    int threads=1;

    TV last_linear_momentum;
    typename TV::SPIN last_angular_momentum;
    T last_te=0;
    T last_grid_ke=0;
    bool output_structures_each_frame=false;
    T quad_F_coeff=0;
    bool asymmetric_system=false;

    VECTOR<BC_TYPE,TV::m*2> side_bc_type; // -x, +x, -y, +y, -z, +z
    std::function<TV(const TV& X,T)> bc_velocity=0;

    MPM_EXAMPLE(const STREAM_TYPE stream_type_input);
    MPM_EXAMPLE(const MPM_EXAMPLE&) = delete;
    void operator=(const MPM_EXAMPLE&) = delete;
    virtual ~MPM_EXAMPLE();
    
    virtual void Write_Output_Files() final;
    virtual void Read_Output_Files() final;
    virtual void Initialize()=0;
    ARRAY<std::function<void(int frame)> > begin_frame;
    ARRAY<std::function<void(int frame)> > end_frame;
    ARRAY<std::function<void()> > write_output_files;
    ARRAY<std::function<void()> > read_output_files;

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
    template<class OBJECT> enable_if_t<!is_pointer<OBJECT>::value>
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
    void Add_Callbacks(bool is_begin,const char* func_name,std::function<void()> func);
//#####################################################################
};
}
#endif
