//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_EXAMPLE__
#define __PLS_EXAMPLE__
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
namespace PhysBAM{

template<class TV> class LEVELSET_MULTIPLE;

template<class TV_input>
class PLS_EXAMPLE:public LEVELSET_CALLBACKS<TV_input>
{
    typedef TV_input TV;
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
    PROJECTION_DYNAMICS_UNIFORM<TV> projection;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV> particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<TV> incompressible;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T> advection_scalar;
    BOUNDARY<TV,T> boundary_scalar;
    BOUNDARY<TV,T> *boundary,*phi_boundary;
    T_BOUNDARY_PHI_WATER phi_boundary_water;
    //ARRAY<T,TV_INT> density,temperature;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV> collision_bodies_affecting_fluid;    

    PLS_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~PLS_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const override
    {V_levelset=face_velocities;}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) override;
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<TV>& levelset_multiple,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const override {}
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Set_Boundary_Conditions(const T time)=0;
    virtual void Adjust_Phi_With_Sources(const T time)=0;
    virtual void Initialize_Phi()=0;
    virtual void Advect_Particles(const T dt,const T time) {}

//#####################################################################
};
}
#endif
