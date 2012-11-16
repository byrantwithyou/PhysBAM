//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPRESSIBLE_EXAMPLE__
#define __COMPRESSIBLE_EXAMPLE__
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Compressible_Fluids/COMPRESSIBLE_FLUID_COLLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
namespace PhysBAM{

template<class TV_input>
class COMPRESSIBLE_EXAMPLE:public EXAMPLE<TV_input>
{
    typedef TV_input TV;
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;typedef VECTOR<T,GRID<TV>::dimension+2> TV_DIMENSION;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename REBIND<T_FACE_ARRAYS_SCALAR,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef EXAMPLE<TV> BASE;
    enum workaround1{d=TV::m};

    using BASE::stream_type;using BASE::initial_time;
    using BASE::output_directory;

public:
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    int number_of_ghost_cells;
    T cfl;

    GRID<TV> mac_grid;
    EULER_UNIFORM<GRID<TV> > euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV> euler_solid_fluid_coupling_utilities;
    COMPRESSIBLE_FLUID_COLLECTION<GRID<TV> > compressible_fluid_collection;
    T_FACE_ARRAYS_SCALAR face_velocities;

    CONSERVATION<GRID<TV>,GRID<TV>::dimension+2>* conservation_method;
    BOUNDARY_REFLECTION_UNIFORM<GRID<TV>,VECTOR<T,TV::m+2> >* boundary;
    BOUNDARY_UNIFORM<GRID<TV>,T>* pressure_boundary;

    T_GRID_BASED_COLLISION_GEOMETRY* collision_bodies_affecting_fluid;

    bool left_wall,right_wall,bottom_wall,top_wall,front_wall,back_wall;
    bool timesplit;
    bool set_max_time_step;T max_time_step;
    int spatial_order;
    int rungekutta_order;
    T tolerance;
    int iterations;
    bool solve_single_cell_neumann_regions;
    bool fluid_affects_solid,solid_affects_fluid;
    bool apply_isobaric_fix;
    bool monitor_conservation_error;
    bool perform_rungekutta_for_implicit_part;
    bool use_sound_speed_for_cfl;
    bool use_sound_speed_based_dt_multiple_for_cfl;
    T multiplication_factor_for_sound_speed_based_dt;

    COMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input);
    virtual ~COMPRESSIBLE_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    virtual void Apply_Isobaric_Fix(const T dt,const T time);
    void Set_Domain_Boundary_Conditions();
    virtual void Set_Dirichlet_Boundary_Conditions(const T time);
    virtual void Set_Neumann_Boundary_Conditions();
    virtual void Set_Boundary_Conditions(const T time);
    virtual void Initialize_Solid_Fluid_Coupling();
    void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
    virtual void Advance_Kinematic_Collision_Bodies(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Initialize_Euler_State()=0;
    virtual void Initialize_Bodies()=0;

//#####################################################################
};
}
#endif
