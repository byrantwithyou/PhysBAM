//#####################################################################
// Copyright 2006-2007, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_SMOKE
//#####################################################################
#ifndef __STANDARD_TESTS_SMOKE__
#define __STANDARD_TESTS_SMOKE__

#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_SMOKE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;typedef UNIFORM_GRID_ITERATOR_FACE<TV> FACE_ITERATOR;typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    using BASE::restart_frame;using BASE::resolution;

    SMOKE_STANDARD_TESTS_3D<GRID<TV> > tests;

    STANDARD_TESTS_SMOKE(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,0,fluids_parameters.SMOKE),
        tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection)
    {
        *fluids_parameters.grid=tests.grid;
    }

    ~STANDARD_TESTS_SMOKE()
    {}

    // Unused callbacks
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}

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
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
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
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    if(tests.test_number!=3||time<tests.explosion_end_time)
        BASE::Adjust_Density_And_Temperature_With_Sources(tests.source,tests.world_to_source,tests.rho,fluids_parameters.temperature_products);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    BASE::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    tests.Initialize_Bodies();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<3> >& face_velocities,ARRAY<bool,FACE_INDEX<3> >& psi_N,const T time) PHYSBAM_OVERRIDE
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
void Get_Divergence(ARRAY<T,VECTOR<int,3> >& divergence,const T dt,const T time) PHYSBAM_OVERRIDE
{
    tests.Get_Divergence(divergence,dt,time);
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const T dt,const T time) PHYSBAM_OVERRIDE
{
    BASE::Get_Body_Force(force,dt,time);
    tests.Get_Body_Force(force,dt,time);
}
//#####################################################################
};
}
#endif
