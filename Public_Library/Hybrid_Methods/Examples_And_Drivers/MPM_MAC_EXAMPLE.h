//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_MAC_EXAMPLE__
#define __MPM_MAC_EXAMPLE__
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/PHASE_ID.h>
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

template<class TV>
class MPM_MAC_EXAMPLE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
    GRID<TV> grid;
    STREAM_TYPE stream_type;

    struct PHASE
    {
        ARRAY<T,FACE_INDEX<TV::m> > mass,volume;
        ARRAY<T,FACE_INDEX<TV::m> > velocity,velocity_save;

        ARRAY<int> valid_flat_indices;
        ARRAY<FACE_INDEX<TV::m> > valid_indices;
        ARRAY<int> simulated_particles;
        GATHER_SCATTER<TV>* gather_scatter;

        PHASE();
        PHASE(const PHASE&) = delete;
        ~PHASE();
        void Initialize(const GRID<TV>& grid,
                        const VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m>& weights,
                        int ghost,int threads);
    };

    ARRAY<PHASE,PHASE_ID> phases;

    // particle stuff
    MPM_PARTICLES<TV>& particles;

    // grid stuff
    ARRAY<TV,FACE_INDEX<TV::m> > location;
    MPM_PROJECTION_SYSTEM<TV>& projection_system;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    MPM_PROJECTION_VECTOR<TV>& sol;
    MPM_PROJECTION_VECTOR<TV>& rhs;
    int ghost;

    // signed distance field & level sets
    ARRAY<T,TV_INT> phi;
    LEVELSET<TV>* levelsets;

    // transfer stuff
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights;
    bool use_affine;
    bool use_early_gradient_transfer;
    T flip;

    // objects
    ARRAY<MPM_COLLISION_OBJECT<TV>*> collision_objects;
    ARRAY<IMPLICIT_OBJECT<TV>* > fluid_walls;

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

    // parameters
    TV gravity;
    T cfl;
    T solver_tolerance;
    int solver_iterations;
    int threads;
    bool use_particle_volumes;
    bool use_preconditioner;

    bool use_shrink,use_reinit;
    T dilation;

    // debugging
    DEBUG_PARTICLES<TV>& debug_particles;
    bool print_stats;
    TV last_linear_momentum;
    typename TV::SPIN last_angular_momentum;
    T last_te;
    T last_grid_ke;
    bool test_system;

    MPM_MAC_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~MPM_MAC_EXAMPLE();
    
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    std::function<void(int frame)> begin_frame;
    std::function<void(int frame)> end_frame;
    std::function<void(T time)> begin_time_step;
    std::function<void(T time)> end_time_step;

    T Potential_Energy(const T time) const;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const;
    void Set_Weights(int order);
    void Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0);
    void Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io);
    template<class OBJECT> typename enable_if<!is_pointer<OBJECT>::value>::type
    Add_Collision_Object(const OBJECT& object,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0)
    {Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<OBJECT>(object),type,friction,func_frame,func_twist);}

    TV Total_Particle_Linear_Momentum() const;
    TV Total_Particle_Linear_Momentum(const PHASE& ph) const;
    TV Total_Grid_Linear_Momentum() const;
    TV Total_Grid_Linear_Momentum(const PHASE& ph) const;
    typename TV::SPIN Total_Grid_Angular_Momentum(T dt) const;
    typename TV::SPIN Total_Grid_Angular_Momentum(const PHASE& ph,T dt) const;
    typename TV::SPIN Total_Particle_Angular_Momentum() const;
    typename TV::SPIN Total_Particle_Angular_Momentum(const PHASE& ph) const;
    T Total_Grid_Kinetic_Energy() const;
    T Total_Grid_Kinetic_Energy(const PHASE& ph) const;
    T Total_Particle_Kinetic_Energy() const;
    T Total_Particle_Kinetic_Energy(const PHASE& ph) const;
    T Average_Particle_Mass() const;
    T Average_Particle_Mass(const PHASE& ph) const;
//#####################################################################
};
}
#endif
