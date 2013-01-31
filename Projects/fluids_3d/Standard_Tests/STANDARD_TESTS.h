//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;typedef UNIFORM_GRID_ITERATOR_FACE<TV> FACE_ITERATOR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;using BASE::Adjust_Phi_With_Source;

    WATER_STANDARD_TESTS_3D<GRID<TV> > tests;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,1,fluids_parameters.WATER),
        tests(*this,fluids_parameters,solid_body_collection.rigid_body_collection)
    {
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    tests.Initialize(test_number,resolution);
    *fluids_parameters.grid=tests.grid;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    tests.Initialize_Advection();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=tests.Initial_Phi(iterator.Location());
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    tests.Initialize_Bodies();
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Update_Fluid_Parameters(dt,time);
    tests.Update_Sources(time);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    for(int s=0;s<tests.sources.m;s++)Adjust_Phi_With_Source(tests.sources(s),tests.world_to_source(s));
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
{
    bool first=true;
    for(int s=0;s<tests.sources.m;s++){Get_Source_Reseed_Mask(tests.sources(s),tests.world_to_source(s),cell_centered_mask,first);first=false;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    for(int s=0;s<tests.sources.m;s++)Get_Source_Velocities(tests.sources(s),tests.world_to_source(s),tests.source_velocity(s));
}
//#####################################################################
// Function Get_Variable_Viscosity
//#####################################################################
void Get_Variable_Viscosity(ARRAY<T,VECTOR<int,3> >& viscosity,const T time) PHYSBAM_OVERRIDE
{
    tests.Get_Variable_Viscosity(viscosity,time);
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    tests.Set_Kinematic_Positions(frame,time,id);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    return tests.Set_Kinematic_Velocities(twist,time,id);
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    tests.Limit_Dt(dt,time);
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() PHYSBAM_OVERRIDE
{
    tests.Initialize_SPH_Particles();
}
//#####################################################################
// Function Get_Analytic_Velocities
//#####################################################################
void Get_Analytic_Velocities(const T time) const PHYSBAM_OVERRIDE
{
    PHYSBAM_FATAL_ERROR("broken");
#if 0
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=tests.Analytic_Velocity(time,iterator.Location())[iterator.Axis()];
#endif
}
//#####################################################################
};
}
#endif
