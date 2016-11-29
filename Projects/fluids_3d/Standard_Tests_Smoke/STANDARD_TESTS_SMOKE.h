//#####################################################################
// Copyright 2006-2007, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_SMOKE
//#####################################################################
#ifndef __STANDARD_TESTS_SMOKE__
#define __STANDARD_TESTS_SMOKE__

#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_SMOKE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::restart_frame;using BASE::resolution;

    SMOKE_STANDARD_TESTS_3D<TV> tests;

    STANDARD_TESTS_SMOKE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.SMOKE),
        tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection)
    {
        parse_args.Parse();
        tests.Initialize(test_number,resolution);
        *fluids_parameters.grid=tests.grid;
    }

    ~STANDARD_TESTS_SMOKE()
    {}

    // Unused callbacks
    void Limit_Dt(T& dt,const T time) override {}
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Preprocess_Frame(const int frame) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Postprocess_Frame(const int frame) override {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
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
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) override
{
    if(tests.test_number!=3||time<tests.explosion_end_time)
        BASE::Adjust_Density_And_Temperature_With_Sources(tests.source,tests.world_to_source,tests.rho,fluids_parameters.temperature_products);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    BASE::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    tests.Initialize_Bodies();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) override
{
    if(tests.Use_Source())
        BASE::Get_Source_Velocities(tests.source,tests.world_to_source,tests.source_velocity);
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    BASE::Construct_Levelsets_For_Objects(time);
}
//#####################################################################
// Function Get_Divergence
//#####################################################################
void Get_Divergence(ARRAY<T,VECTOR<int,3> >& divergence,const T dt,const T time) override
{
    tests.Get_Divergence(divergence,dt,time);
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time) override
{
    BASE::Get_Body_Force(force,dt,time);
    tests.Get_Body_Force(force,dt,time);
}
//#####################################################################
};
}
#endif
