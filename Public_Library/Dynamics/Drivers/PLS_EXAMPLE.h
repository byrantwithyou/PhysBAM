//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_EXAMPLE__
#define __PLS_EXAMPLE__
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_BODY_COLLECTION_POLICY_UNIFORM.h>
#include <Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class LEVELSET_MULTIPLE;

template<class TV_input>
class PLS_EXAMPLE:public LEVELSET_CALLBACKS<GRID<TV_input> >
{
    typedef TV_input TV;
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef BOUNDARY_PHI_WATER<TV> T_BOUNDARY_PHI_WATER;
    typedef typename COLLISION_BODY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
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
    MPI_UNIFORM_GRID<GRID<TV> > *mpi_grid;
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > projection;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> > particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<GRID<TV> > incompressible;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T> advection_scalar;
    BOUNDARY<TV,T> boundary_scalar;
    BOUNDARY<TV,T> *boundary,*phi_boundary;
    T_BOUNDARY_PHI_WATER phi_boundary_water;
    //ARRAY<T,TV_INT> density,temperature;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    T_GRID_BASED_COLLISION_GEOMETRY collision_bodies_affecting_fluid;    

    PLS_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~PLS_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE
    {V_levelset=face_velocities;}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE
    {
        if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

        TV& X=particles.X(index);TV X_new=X+dt*V;
        T max_collision_distance=particle_levelset_evolution.Particle_Levelset(0).Particle_Collision_Distance(particles.quantized_collision_distance(index));
        T min_collision_distance=particle_levelset_evolution.Particle_Levelset(0).min_collision_distance_factor*max_collision_distance;
        TV min_corner=mac_grid.domain.Minimum_Corner(),max_corner=mac_grid.domain.Maximum_Corner();
        for(int axis=0;axis<GRID<TV>::dimension;axis++){
            if(domain_boundary[axis][0] && X_new[axis]<min_corner[axis]+max_collision_distance){
                T collision_distance=X[axis]-min_corner[axis];
                if(collision_distance>max_collision_distance)collision_distance=X_new[axis]-min_corner[axis];
                collision_distance=max(min_collision_distance,collision_distance);
                X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
                V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
            if(domain_boundary[axis][1] && X_new[axis]>max_corner[axis]-max_collision_distance){
                T collision_distance=max_corner[axis]-X[axis];
                if(collision_distance>max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
                collision_distance=max(min_collision_distance,collision_distance);
                X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
                V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
    }
    
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<GRID<TV> >& levelset_multiple,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE {}
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
