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
#include <functional>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
template<class TV> class GATHER_SCATTER;
template<class TV> class IMPLICIT_OBJECT;
template<class TV> class MPM_PARTICLES;
template<class TV> class PARTICLE_GRID_FORCES;
template<class TV> class PARTICLE_GRID_WEIGHTS;
template<class TV> class MPM_PLASTICITY_MODEL;
template<class TV> class PROJECTION_UNIFORM;

template<class TV>
class MPM_MAC_EXAMPLE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
    GRID<TV> grid;
    STREAM_TYPE stream_type;

    // particle stuff
    MPM_PARTICLES<TV>& particles;
    ARRAY<int> simulated_particles;
    ARRAY<bool> particle_is_simulated;

    // grid stuff
    ARRAY<T,FACE_INDEX<TV::m> > mass,volume,density;
    ARRAY<T,FACE_INDEX<TV::m> > velocity;
    ARRAY<TV,FACE_INDEX<TV::m> > location;
    ARRAY<int> valid_flat_indices;
    ARRAY<FACE_INDEX<TV::m> > valid_indices;
    PROJECTION_UNIFORM<TV>* projection;
    int ghost;

    // transfer stuff
    VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights;
    GATHER_SCATTER<TV>& gather_scatter;
    bool use_affine;
    bool use_early_gradient_transfer;

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
    int restart;
    T dt,time,frame_dt,min_dt,max_dt;
    bool only_write_particles;

    // parameters
    TV gravity;
    T cfl;
    T solver_tolerance;
    int solver_iterations;
    int threads;

    // debugging
    DEBUG_PARTICLES<TV>& debug_particles;
    bool print_stats;
    TV last_linear_momentum;
    typename TV::SPIN last_angular_momentum;
    T last_te;
    T last_grid_ke;

    MPM_MAC_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~MPM_MAC_EXAMPLE();
    
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;
    virtual void Begin_Frame(const int frame)=0;
    virtual void End_Frame(const int frame)=0;
    virtual void Begin_Time_Step(const T time)=0;
    virtual void End_Time_Step(const T time)=0;

    T Potential_Energy(const T time) const;
    void Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const;
    void Set_Weights(int order,int threads);
    void Add_Collision_Object(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0);
    void Add_Fluid_Wall(IMPLICIT_OBJECT<TV>* io);
    template<class OBJECT> typename enable_if<!is_pointer<OBJECT>::value>::type
    Add_Collision_Object(const OBJECT& object,COLLISION_TYPE type,T friction,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0)
    {Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<OBJECT>(object),type,friction,func_frame,func_twist);}

    TV Total_Particle_Linear_Momentum() const;
    TV Total_Grid_Linear_Momentum() const;
    typename TV::SPIN Total_Grid_Angular_Momentum(T dt) const;
    typename TV::SPIN Total_Particle_Angular_Momentum() const;
    T Total_Grid_Kinetic_Energy() const;
    T Total_Particle_Kinetic_Energy() const;
    T Average_Particle_Mass() const;
//#####################################################################
};
}
#endif
