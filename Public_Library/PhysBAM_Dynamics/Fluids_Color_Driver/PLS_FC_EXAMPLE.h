//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_FC_EXAMPLE__
#define __PLS_FC_EXAMPLE__
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_EXTRAPOLATE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Finite_Elements/LEVELSET_COLOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class LEVELSET_MULTIPLE_UNIFORM;

template<class TV>
class PLS_FC_EXAMPLE:public LEVELSET_CALLBACKS<GRID<TV> >
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename LEVELSET_POLICY<GRID<TV> >::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
    typedef BOUNDARY_PHI_WATER<GRID<TV> > T_BOUNDARY_PHI_WATER;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    T initial_time;
    int last_frame;
    std::string frame_title;
    int write_substeps_level;
    bool write_output_files;
    std::string output_directory;
    int restart;
    int number_of_ghost_cells;

    T dt;
    int time_steps_per_frame;
    bool use_preconditioner;
    int max_iter;
    bool dump_matrix;
    bool wrap;
    ARRAY<T> mu;

    GRID<TV> grid;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> > particle_levelset_evolution;
    ARRAY<int,FACE_INDEX<TV::dimension> > face_color;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T> advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary,*phi_boundary;
    BOUNDARY_EXTRAPOLATE_CELL<TV,T> cell_extrapolate;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    LEVELSET_COLOR<TV> levelset_color;
    T_GRID_BASED_COLLISION_GEOMETRY collision_bodies_affecting_fluid;

    PLS_FC_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~PLS_FC_EXAMPLE();
    
    void Get_Levelset_Velocity(const GRID<TV>& grid,T_LEVELSET& levelset,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE
    {V_levelset=face_velocities;}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE_UNIFORM<GRID<TV> >& levelset_multiple,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE {}
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Initialize()=0;

//#####################################################################
};
}
#endif
