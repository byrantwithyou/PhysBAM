//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPRESSIBLE_EXAMPLE__
#define __COMPRESSIBLE_EXAMPLE__
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Compressible/Compressible_Fluids/COMPRESSIBLE_FLUID_COLLECTION.h>
#include <Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
namespace PhysBAM{

template<class TV_input>
class COMPRESSIBLE_EXAMPLE:public EXAMPLE<TV_input>
{
    typedef TV_input TV;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef EXAMPLE<TV> BASE;
    enum workaround1{d=TV::m};

    using BASE::stream_type;

public:
    using BASE::last_frame;using BASE::frame_rate;using BASE::viewer_dir;
    int number_of_ghost_cells;
    T cfl;

    GRID<TV> mac_grid;
    EULER_UNIFORM<TV> euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV> euler_solid_fluid_coupling_utilities;
    COMPRESSIBLE_FLUID_COLLECTION<TV> compressible_fluid_collection;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities;

    CONSERVATION<TV,TV::m+2>* conservation_method;
    BOUNDARY_REFLECTION_UNIFORM<TV,VECTOR<T,TV::m+2> >* boundary;
    BOUNDARY<TV,T>* pressure_boundary;

    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_bodies_affecting_fluid;

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

    COMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~COMPRESSIBLE_EXAMPLE();
    
    T Time_At_Frame(const int frame) const override
    {return frame/frame_rate;}

    virtual void Apply_Isobaric_Fix(const T dt,const T time);
    void Set_Domain_Boundary_Conditions();
    virtual void Set_Dirichlet_Boundary_Conditions(const T time);
    virtual void Set_Neumann_Boundary_Conditions();
    virtual void Set_Boundary_Conditions(const T time);
    virtual void Initialize_Solid_Fluid_Coupling();
    void Write_Output_Files() const override;
    virtual void Advance_Kinematic_Collision_Bodies(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Initialize_Euler_State()=0;
    virtual void Initialize_Bodies()=0;

//#####################################################################
};
}
#endif
