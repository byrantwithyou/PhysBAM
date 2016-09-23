//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_MULTIPHASE
//#####################################################################
#ifndef __STANDARD_TESTS_MULTIPHASE__
#define __STANDARD_TESTS_MULTIPHASE__

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_MULTIPHASE_2D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_MULTIPHASE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;using BASE::fluid_collection;
    using BASE::solid_body_collection;using BASE::test_number;using BASE::resolution;using BASE::Adjust_Phi_With_Source;

    WATER_STANDARD_TESTS_MULTIPHASE_2D<TV> tests;

    STANDARD_TESTS_MULTIPHASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.WATER),
        tests(*this,fluids_parameters,fluid_collection,solid_body_collection.rigid_body_collection)
    {
        parse_args.Parse();
        fluids_parameters.Initialize_Number_Of_Regions(WATER_STANDARD_TESTS_MULTIPHASE_2D<TV>::Number_Of_Regions(test_number));
        tests.Initialize(test_number,resolution);
        *fluids_parameters.grid=tests.grid;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.solve_neumann_regions=true;
    }

    ~STANDARD_TESTS_MULTIPHASE()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) override {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Postprocess_Frame(const int frame) override {}
    void Postprocess_Phi(const T time) override {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initial_Phi
//#####################################################################
void Initialize_Advection() override
{
    tests.Initialize_Advection();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    for(int i=0;i<fluids_parameters.number_of_regions;i++)for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluids_parameters.particle_levelset_evolution_multiple->phis(i)(iterator.Cell_Index())=tests.Initial_Phi(i,iterator.Location());        
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Velocities() override
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(iterator.Face_Index())=tests.Initial_Velocity(iterator.Location())[axis];}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    tests.Initialize_Bodies();
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
    for(int s=0;s<tests.sources.m;s++)Adjust_Phi_With_Source(tests.sources(s),tests.source_region(s),tests.world_to_source(s));
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
    if(tests.test_number==16&&fabs(fluids_parameters.grid->domain.min_corner.y)<1e-5){
        ARRAY_VIEW<T,VECTOR<int,2> >& u=fluid_collection.incompressible_fluid_collection.face_velocities.Component(1);
        ARRAY_VIEW<bool,VECTOR<int,2> >& psi_N_u=fluids_parameters.incompressible_multiphase->projection.elliptic_solver->psi_N.Component(1);
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid,1,GRID<TV>::GHOST_REGION,2);iterator.Valid();iterator.Next()){
            //if(fluids_parameters.particle_levelset_evolution->Levelset(1).phi(iterator.Cell_Index())<=0){
            VECTOR<int,2> cell=iterator.Cell_Index();
            if(tests.armadillo->phi(cell+VECTOR<int,2>(0,1))<=0){
                if(u.Valid_Index(cell+VECTOR<int,2>(0,1)))u(cell+VECTOR<int,2>(0,1))=0;
                if(u.Valid_Index(cell+VECTOR<int,2>(1,1)))u(cell+VECTOR<int,2>(1,1))=0;
                if(psi_N_u.Valid_Index(cell+VECTOR<int,2>(0,1)))psi_N_u(cell+VECTOR<int,2>(0,1))=true;
                if(psi_N_u.Valid_Index(cell+VECTOR<int,2>(1,1)))psi_N_u(cell+VECTOR<int,2>(1,1))=true;}}}

    for(int s=0;s<tests.sources.m;s++)Get_Source_Velocities(tests.sources(s),tests.world_to_source(s),tests.source_velocity(s));
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) override
{
    tests.Limit_Dt(dt,time);
}
//#####################################################################
};
}
#endif
