//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MAC_EXAMPLE__
#define __MPM_MAC_EXAMPLE__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
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
    GRID<TV> grid;
    STREAM_TYPE stream_type;
    enum BC_TYPE {BC_FREE, BC_SLIP, BC_NOSLIP, BC_PERIODIC};
    
    ARRAY<T,FACE_INDEX<TV::m> > mass,volume;
    ARRAY<T,FACE_INDEX<TV::m> > velocity,velocity_save;

    ARRAY<int> valid_flat_indices;
    ARRAY<FACE_INDEX<TV::m> > valid_indices;
    ARRAY<int> simulated_particles;
    GATHER_SCATTER<TV>* gather_scatter=0;
    T density=0; // if not using per-particle mass
    T viscosity=0;

    // signed distance field & level sets
    ARRAY<T,TV_INT> phi;
    LEVELSET<TV>* levelset=0;

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
    int ghost;
    VECTOR<BC_TYPE,TV::m*2> bc_type; // -x, +x, -y, +y, -z, +z
    // Valid if BC_SLIP or BC_NOSLIP; velocity at face. null=0
    VECTOR<std::function<T(const TV& X,int axis,T)>,TV::m*2> bc_velocity;
    // Valid if BC_FREE; pressure at ghost cell. null=0
    std::function<T(TV_INT,T)> bc_pressure;
    BOUNDARY_MAC_GRID_PERIODIC<TV,T>& periodic_boundary;
    ARRAY<int,TV_INT> cell_index;
    int dof;
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
    char extrap_type;
    bool clamp_particles;

    // transfer stuff
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights;
    bool use_affine;
    bool lag_Dp;
    T flip;

    // objects
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    ARRAY<IMPLICIT_OBJECT<TV>* > fluid_walls;
    RANDOM_NUMBERS<T> random;

    // initialization & output
    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    int substeps_delay_frame;
    std::string output_directory;
    std::string data_directory;
    std::string test_output_prefix;
    bool use_test_output;
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    bool only_write_particles;
    bool only_log;

    // parameters
    TV gravity;
    T cfl;
    T solver_tolerance;
    int solver_iterations;
    int threads;
    bool use_particle_volumes;
    bool use_constant_density;
    bool move_mass_inside;
    bool move_mass_inside_nearest;
    bool use_preconditioner;
    bool use_phi;
    int rk_particle_order;
    bool use_massless_particles;
    bool use_multiphase_projection;
    bool use_bump;
    bool use_reseeding;
    T radius_sphere,radius_escape;
    bool use_periodic_test_shift;
    TV_INT periodic_test_shift;
    // d: default; c: always use pic update; p: use pic in the first step only
    char position_update;
    
    // debugging
    DEBUG_PARTICLES<TV>& debug_particles;
    TV last_linear_momentum;
    typename TV::SPIN last_angular_momentum;
    T last_grid_te=0;
    T last_grid_ke=0;
    T last_part_te=0;
    T last_part_ke=0;
    bool test_system;
    bool print_matrix;
    bool particle_vort;
    bool use_object_extrap=false;

    MPM_MAC_EXAMPLE(const STREAM_TYPE stream_type_input);
    MPM_MAC_EXAMPLE(const MPM_MAC_EXAMPLE&) = delete;
    void operator=(const MPM_MAC_EXAMPLE&) = delete;
    virtual ~MPM_MAC_EXAMPLE();
    
    void Write_Output_Files(const int frame);
    void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    ARRAY<std::function<void(int frame)> > begin_frame;
    ARRAY<std::function<void(int frame)> > end_frame;
    ARRAY<std::function<void(int frame)> > write_output_files;
    ARRAY<std::function<void(int frame)> > read_output_files;

    HASHTABLE<std::string,PAIR<bool,VECTOR<ARRAY<std::function<void()> >,2> > > time_step_callbacks; // begin, end

    T Potential_Energy(const T time) const;
    T Total_Particle_Vorticity() const;
    void Apply_Forces(const T time);
    virtual TV Compute_Analytic_Force(const TV& X,const T time) const;
    void Set_Weights(int order);
    void Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0);
    void Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io);
    template<class OBJECT> typename enable_if<!is_pointer<OBJECT>::value>::type
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
