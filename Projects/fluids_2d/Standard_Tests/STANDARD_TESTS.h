//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,2> TV;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::test_number;using BASE::fluid_collection;using BASE::solid_body_collection;using BASE::resolution;using BASE::Adjust_Phi_With_Source;

    WATER_STANDARD_TESTS_2D<TV> tests;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,1,fluids_parameters.WATER),
        tests(*this,fluids_parameters,solid_body_collection.rigid_body_collection)
    {
        parse_args.Parse();
        tests.Initialize(test_number,resolution);
        *fluids_parameters.grid=tests.grid;
        fluids_parameters.write_ghost_values=true;
        fluids_parameters.store_particle_ids=true;
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) override {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) override {}
    void Postprocess_Phi(const T time) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Extrapolate_Phi_Into_Objects(const T time) override {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    tests.Initialize_Advection();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() override
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=tests.Initial_Phi(iterator.Location());
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    tests.Initialize_Bodies();
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    tests.Update_Rigid_Bodies(time);
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Construct_Levelsets_For_Objects(time);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
    tests.Update_Sources(time);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    for(int s=0;s<tests.sources.m;s++)Adjust_Phi_With_Source(tests.sources(s),tests.world_to_source(s));
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time) override
{
    bool first=true;
    for(int s=0;s<tests.sources.m;s++){Get_Source_Reseed_Mask(tests.sources(s),tests.world_to_source(s),cell_centered_mask,first);first=false;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override
{
    for(int s=0;s<tests.sources.m;s++)Get_Source_Velocities(tests.sources(s),tests.world_to_source(s),tests.source_velocity(s));
}
//#####################################################################
// Function Get_Variable_Viscosity
//#####################################################################
void Get_Variable_Viscosity(ARRAY<T,VECTOR<int,2> >& viscosity,const T time) override
{
    tests.Get_Variable_Viscosity(viscosity,time);
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    tests.Set_Kinematic_Positions(frame,time,id);
    BASE::Set_Kinematic_Positions(frame,time,id);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override
{
    return BASE::Set_Kinematic_Velocities(twist,time,id) || tests.Set_Kinematic_Velocities(twist,time,id);
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) override
{
    tests.Limit_Dt(dt,time);
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() override
{
    tests.Initialize_SPH_Particles();
}
//#####################################################################
// Function Get_Analytic_Velocities
//#####################################################################
void Get_Analytic_Velocities(const T time) const override {
    PHYSBAM_FATAL_ERROR("broken");
#if 0
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=tests.Analytic_Velocity(time,iterator.Location())[iterator.Axis()];
#endif
}
//#####################################################################
};
}
#endif
