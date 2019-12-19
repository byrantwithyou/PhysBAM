//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MAC_EXAMPLE__
#define __MPM_MAC_EXAMPLE__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <functional>
namespace PhysBAM{

template<class TV> class LEVELSET;
template<class TV> class DEBUG_PARTICLES;
template<class TV> class GATHER_SCATTER;
template<class TV> class IMPLICIT_OBJECT;
template<class TV> class MPM_PARTICLES;
template<class TV> class PARTICLE_GRID_FORCES;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_PLASTICITY_MODEL;
template<class TV> class PROJECTION_UNIFORM;
template<class TV> class MPM_PROJECTION_SYSTEM;
template<class TV> class MPM_PROJECTION_VECTOR;
template<class TV> class KRYLOV_VECTOR_BASE;
template<class TV,class T2> class BOUNDARY_MAC_GRID_PERIODIC;

template<class TV>
class MPM_MAC_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;

    STREAM_TYPE stream_type;
    int ghost=3;
    GRID<TV> grid;
    enum BC_TYPE {BC_FREE, BC_SLIP, BC_NOSLIP, BC_PERIODIC};
    
    ARRAY<T,FACE_INDEX<TV::m> > mass,mass_save,volume;
    ARRAY<T,FACE_INDEX<TV::m> > velocity,velocity_save;
    ARRAY<bool,FACE_INDEX<TV::m> > valid_xfer_data;

    ARRAY<T,TV_INT> pressure_save;
    ARRAY<bool,TV_INT> pressure_valid;
    bool use_warm_start=false;
    
    ARRAY<int> valid_flat_indices;
    ARRAY<FACE_INDEX<TV::m> > valid_indices;
    ARRAY<int> simulated_particles;
    GATHER_SCATTER<TV>* gather_scatter=0;
    T density=0; // if not using per-particle mass
    T viscosity=0;

    // signed distance field & level sets
    GRID<TV> levelset_grid;
    ARRAY<T,TV_INT> phi;
    LEVELSET<TV>& levelset;
    bool disable_free_surface=false;

    // XPIC stuff
    int xpic=0;
    ARRAY<T,FACE_INDEX<TV::m> > xpic_v,xpic_v_star;
    ARRAY_VIEW<TV> effective_v;

    // particle stuff
    MPM_PARTICLES<TV>& particles;
    VECTOR<ARRAY<SYMMETRIC_MATRIX<T,TV::m> >,TV::m> Dp_inv;

    // grid stuff
    ARRAY<TV,FACE_INDEX<TV::m> > location;
    MPM_PROJECTION_SYSTEM<TV>& projection_system;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    MPM_PROJECTION_VECTOR<TV>& sol;
    MPM_PROJECTION_VECTOR<TV>& rhs;
    VECTOR<BC_TYPE,TV::m*2> side_bc_type; // -x, +x, -y, +y, -z, +z
    std::function<BC_TYPE(const TV& X,T)> bc_type=0;
    // Valid if BC_SLIP or BC_NOSLIP; velocity at face. null=0
    std::function<TV(const TV& X,T)> bc_velocity=0;
    // Valid if BC_FREE; pressure at ghost cell. null=0
    std::function<T(TV,T)> bc_pressure=0;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T>& periodic_boundary;
    ARRAY<int,TV_INT> cell_index;
    int dof=0;
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N;
    ARRAY<T,FACE_INDEX<TV::m> > force;
    // Extrapolation functions.
    // Options:
    // 0: fill zeros
    // C: copy nearest interior
    // c: copy nearest interior (use BC values if any)
    // L: linearly extrapolate from interior
    // l: linearly extrapolate from interior (use BC values if any)
    // a: extrpolate using analytic function
    // r: reflect the interior for tangential directions; reflect with negation for normal directions
    char extrap_type='p';
    bool clamp_particles=false;

    // transfer stuff
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights;
    bool use_affine=true;
    T flip=0;
    ARRAY_VIEW<TV> flip_adv_velocity;

    // objects
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    ARRAY<IMPLICIT_OBJECT<TV>* > fluid_walls;
    RANDOM_NUMBERS<T> random;

    // initialization & output
    int last_frame=100;
    std::string frame_title;
    int write_substeps_level=-1;
    int substeps_delay_frame=-1;
    VIEWER_DIR viewer_dir{"output"};
    std::string data_directory="../../Public_Data";
    std::string test_output_prefix;
    bool use_test_output=false;
    bool auto_restart=false;
    int restart=0;
    T dt=0,time=0,frame_dt=(T)1/24,min_dt=0,max_dt=(T)1/24;
    bool only_write_particles=false;
    bool only_log=false;

    // parameters
    TV gravity;
    T cfl=1;
    T solver_tolerance=std::numeric_limits<T>::epsilon()*10;
    int solver_iterations=1000;
    int threads=1;
    bool use_particle_volumes=false;
    bool use_constant_density=true;
    bool use_preconditioner=true;
    bool use_phi=false;
    int rk_particle_order=0;
    bool use_massless_particles=false;
    bool use_level_set_projection=false;
    bool use_reseeding=false;
    bool use_periodic_test_shift=false;
    TV_INT periodic_test_shift;
    // d: default; c: always use pic update; p: use pic in the first step only
    char position_update='d';
    
    // debugging
    DEBUG_PARTICLES<TV>& debug_particles;
    TV last_linear_momentum;
    typename TV::SPIN last_angular_momentum;
    T last_grid_te=0;
    T last_grid_ke=0;
    T last_part_te=0;
    T last_part_ke=0;
    bool test_system=false;
    bool print_matrix=false;
    bool particle_vort=false;
    bool use_object_extrap=false;
    bool use_mls_xfers=false;
    bool zero_invalid=false;
    
    MPM_MAC_EXAMPLE(const STREAM_TYPE stream_type_input);
    MPM_MAC_EXAMPLE(const MPM_MAC_EXAMPLE&) = delete;
    void operator=(const MPM_MAC_EXAMPLE&) = delete;
    virtual ~MPM_MAC_EXAMPLE();
    
    void Write_Output_Files();
    void Read_Output_Files();
    virtual void Initialize()=0;
    ARRAY<std::function<void(int frame)> > begin_frame;
    ARRAY<std::function<void(int frame)> > end_frame;
    ARRAY<std::function<void()> > write_output_files;
    ARRAY<std::function<void()> > read_output_files;

    HASHTABLE<std::string,PAIR<bool,VECTOR<ARRAY<std::function<void()> >,2> > > time_step_callbacks; // begin, end

    T Potential_Energy(const T time) const;
    T Total_Particle_Vorticity() const;
    void Apply_Forces(const T time);
    virtual TV Compute_Analytic_Force(const TV& X,const T time) const;
    void Set_Weights(int order);
    void Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0);
    void Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io);
    template<class OBJECT> enable_if_t<!is_pointer<OBJECT>::value>
    Add_Collision_Object(const OBJECT& object,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0)
    {Add_Collision_Object(Make_IO(object),type,friction,func_frame,func_twist);}

    TV Total_Particle_Linear_Momentum() const;
    TV Total_Grid_Linear_Momentum() const;
    typename TV::SPIN Total_Grid_Angular_Momentum(T dt) const;
    typename TV::SPIN Total_Particle_Angular_Momentum() const;
    T Total_Grid_Kinetic_Energy() const;
    T Total_Particle_Kinetic_Energy() const;
    T Average_Particle_Mass() const;
    void Execute_Callbacks(bool is_begin,const char* func_name);
    void Add_Callbacks(bool is_begin,const char* func_name,std::function<void()> func);
    void Print_Grid_Stats(const char* str);
    void Print_Particle_Stats(const char* str);
    void Dump_Grid_ShiftTest(const std::string& var_name,const ARRAY<T,FACE_INDEX<TV::m> >& arr);
//#####################################################################
};
}
#endif
