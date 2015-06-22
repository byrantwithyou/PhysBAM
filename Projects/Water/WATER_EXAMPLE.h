//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WATER_EXAMPLE__
#define __WATER_EXAMPLE__
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Incompressible/Projection/PROJECTION_UNIFORM.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
namespace PhysBAM{

template<class TV> class LEVELSET_MULTIPLE_UNIFORM;

//TODO: Should adventually derive off of a incompressible project
template<class TV>
class WATER_EXAMPLE:public LEVELSET_CALLBACKS<TV>,public RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef BOUNDARY_PHI_WATER<TV> T_BOUNDARY_PHI_WATER;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    std::string frame_title;
    int write_substeps_level;
    bool write_output_files;
    std::string output_directory;
    int restart;
    int number_of_ghost_cells;

    T cfl;

    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<TV> *mpi_grid;
    THREAD_QUEUE* thread_queue;
    PROJECTION_DYNAMICS_UNIFORM<TV>& projection;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV> particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<TV> incompressible;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T> advection_scalar;
    BOUNDARY<TV,T> boundary_scalar;
    BOUNDARY<TV,T> *boundary,*phi_boundary;
    T_BOUNDARY_PHI_WATER phi_boundary_water;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    RIGID_BODY_COLLECTION<TV> rigid_body_collection;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV> collision_bodies_affecting_fluid;
    ARRAY<IMPLICIT_OBJECT<TV>*> sources;

    WATER_EXAMPLE(const STREAM_TYPE stream_type_input,int number_of_threads=1);
    virtual ~WATER_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE
    {V_levelset=face_velocities;}

    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
    //{T range=0.7;int int_time=(int)time;frame.t=TV::All_Ones_Vector()*0.5;
    //if(int_time%2==0) frame.t.x=(time-int_time)*range+(1-range)/2;
    //else frame.t.x=1-(time-int_time)*range-(1-range)/2;}
    {T range=0.6;frame.t=TV::All_Ones_Vector()*0.5;
    if(time<=2) frame.t(2)=time*range+(1-range)/2.;}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
    {return false;}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<TV>& levelset_multiple,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE {}
    void Initialize_Grid(TV_INT counts,RANGE<TV> range);
    void Write_Output_Files(const int frame);
    void Read_Output_Files(const int frame);
    void Set_Boundary_Conditions(const T time);
    void Adjust_Phi_With_Sources(const T time);
    void Adjust_Phi_With_Objects(const T time);
    void Extrapolate_Phi_Into_Objects(const T time);
    void Initialize_Phi();

//#####################################################################
};
}
#endif
