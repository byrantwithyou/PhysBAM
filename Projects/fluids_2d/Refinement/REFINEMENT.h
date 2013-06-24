//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REFINEMENT
//#####################################################################
#ifndef __REFINEMENT__
#define __REFINEMENT__

#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
namespace PhysBAM{

template<class T_input>
class REFINEMENT:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::data_directory;using BASE::fluid_collection;using BASE::solid_body_collection;
    using BASE::resolution;using BASE::test_number;using BASE::parse_args;using BASE::Get_Object_Velocities; // silence -Woverloaded-virtual

    SMOKE_STANDARD_TESTS_2D<GRID<TV> > tests;
    T angle_fraction;

    REFINEMENT(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,0,fluids_parameters.SMOKE),
        tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection),angle_fraction(0)
    {
    }

    virtual ~REFINEMENT()
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
    parse_args->Add("-angle_fraction",&angle_fraction,"fraction","Angle fraction");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    tests.Initialize(test_number,resolution,angle_fraction);
    output_directory="Refinement/output";
    *fluids_parameters.grid=tests.grid;
    last_frame=100;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {}
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
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=tests.Initial_Velocity(iterator.Location())[iterator.Axis()];
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
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
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<2> >& face_velocities,ARRAY<bool,FACE_INDEX<2> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
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
void Get_Divergence(ARRAY<T,VECTOR<int,2> >& divergence,const T dt,const T time) PHYSBAM_OVERRIDE
{
    tests.Get_Divergence(divergence,dt,time);
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE
{
    tests.Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& projection,const T dt,const T time) PHYSBAM_OVERRIDE
{
    BASE::Get_Object_Velocities(projection,dt,time);
    tests.Get_Object_Velocities(projection,dt,time);
}
//#####################################################################
};
}
#endif
