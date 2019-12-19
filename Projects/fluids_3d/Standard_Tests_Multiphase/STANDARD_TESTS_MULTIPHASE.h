//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_MULTIPHASE
//#####################################################################
#ifndef __STANDARD_TESTS_MULTIPHASE__
#define __STANDARD_TESTS_MULTIPHASE__

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_MULTIPHASE_3D.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS_MULTIPHASE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::viewer_dir;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::solid_body_collection;using BASE::Time_At_Frame;using BASE::test_number;using BASE::resolution;using BASE::Adjust_Phi_With_Source;

    WATER_STANDARD_TESTS_MULTIPHASE_3D<TV> tests;

    STANDARD_TESTS_MULTIPHASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.WATER),
        tests(*this,fluids_parameters,fluid_collection,solid_body_collection.rigid_body_collection)
    {
        parse_args.Parse();

        tests.Initialize(test_number,resolution,restart);
        fluids_parameters.Initialize_Number_Of_Regions(WATER_STANDARD_TESTS_MULTIPHASE_3D<TV>::Number_Of_Regions(test_number));
        *fluids_parameters.grid=tests.grid;
/*
  fluids_parameters.viscosities(0)=(T)100;
  fluids_parameters.implicit_viscosity=true;
  fluids_parameters.incompressible_iterations=200;
  fluids_parameters.implicit_viscosity_iterations=200;

  fluids_parameters.use_multiphase_strain(0)=true;
  fluids_parameters.elastic_moduli(0)=10000;
  fluids_parameters.plasticity_alphas(0)=1;
  fluids_parameters.cfl/=4;*/
        this->limit_dt=[this](T& dt,T time)
            {
                T cfl_number=(T).8;
                for(int i=0;i<fluids_parameters.incompressible_multiphase->strains.m;i++)
                    if(fluids_parameters.incompressible_multiphase->strains(i))
                        dt=min(cfl_number*fluids_parameters.incompressible_multiphase->strains(i)->CFL(fluids_parameters.densities(i)),dt);
                tests.Limit_Dt(dt,time);
            };
    }

    ~STANDARD_TESTS_MULTIPHASE()
    {}

//#####################################################################
// Function Initial_Advection
//#####################################################################
void Initialize_Advection()    override
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
    for(int s=0;s<tests.sources.m;s++)Adjust_Phi_With_Source(tests.sources(s),tests.source_region(s),tests.world_to_source(s));

    if(tests.test_number==16&&time>=Time_At_Frame(501)&&time<Time_At_Frame(502)){
        static bool has_added_ball=false;
        if(!has_added_ball){
            Adjust_Phi_With_Source(SPHERE<TV>(VECTOR<T,3>((T).5,(T).8,(T).5),(T).15),4,MATRIX<T,4>::Identity_Matrix());
            has_added_ball=true;}}
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) override
{
    bool first=true;
    for(int s=0;s<tests.sources.m;s++){Get_Source_Reseed_Mask(tests.sources(s),tests.world_to_source(s),cell_centered_mask,first);first=false;}
    if(tests.test_number==16&&time>=Time_At_Frame(501)&&time<Time_At_Frame(502)){
        static bool has_added_ball=false;
        if(!has_added_ball){
            Get_Source_Reseed_Mask(SPHERE<TV>(VECTOR<T,3>((T).5,(T).8,(T).5),(T).15),MATRIX<T,4>::Identity_Matrix(),cell_centered_mask,first);
            has_added_ball=true;}}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override
{
#if 0
    if(tests.test_number==16&&fabs(fluids_parameters.grid->domain.min_corner.y)<1e-5&&time<1.1333333333+1e-5){
        ARRAY<T,VECTOR<int,3> >& u=fluids_parameters.incompressible_multiphase->projection.face_velocities.Component(0);
        ARRAY<T,VECTOR<int,3> >& w=fluids_parameters.incompressible_multiphase->projection.face_velocities.Component(2);
        ARRAY<bool,VECTOR<int,3> >& psi_N_u=fluids_parameters.incompressible_multiphase->projection.elliptic_solver->psi_N.Component(0);
        ARRAY<bool,VECTOR<int,3> >& psi_N_w=fluids_parameters.incompressible_multiphase->projection.elliptic_solver->psi_N.Component(2);
        for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid,1,GRID<TV>::GHOST_REGION,2);iterator.Valid();iterator.Next()){
            //if(fluids_parameters.particle_levelset_evolution->Levelset(0).phi(iterator.Cell_Index())<=0){
            VECTOR<int,3> cell=iterator.Cell_Index();
            if(tests.armadillo->phi(cell+VECTOR<int,3>(0,1,0))<=0){
                if(u.Valid_Index(cell+VECTOR<int,3>(0,1,0)))u(cell+VECTOR<int,3>(0,1,0))=0;
                if(u.Valid_Index(cell+VECTOR<int,3>(1,1,0)))u(cell+VECTOR<int,3>(1,1,0))=0;
                if(w.Valid_Index(cell+VECTOR<int,3>(0,1,0)))w(cell+VECTOR<int,3>(0,1,0))=0;
                if(w.Valid_Index(cell+VECTOR<int,3>(0,1,1)))w(cell+VECTOR<int,3>(0,1,1))=0;
                if(psi_N_u.Valid_Index(cell+VECTOR<int,3>(0,1,0)))psi_N_u(cell+VECTOR<int,3>(0,1,0))=true;
                if(psi_N_u.Valid_Index(cell+VECTOR<int,3>(1,1,0)))psi_N_u(cell+VECTOR<int,3>(1,1,0))=true;
                if(psi_N_w.Valid_Index(cell+VECTOR<int,3>(0,1,0)))psi_N_w(cell+VECTOR<int,3>(0,1,0))=true;
                if(psi_N_w.Valid_Index(cell+VECTOR<int,3>(0,1,1)))psi_N_w(cell+VECTOR<int,3>(0,1,1))=true;}}}

    for(int s=0;s<tests.sources.m;s++)Get_Source_Velocities(tests.sources(s),tests.world_to_source(s),tests.source_velocity(s));

    if(tests.test_number==16&&time>=Time_At_Frame(501)&&time<Time_At_Frame(502)){
        static bool has_added_ball=false;
        if(!has_added_ball){
            Get_Source_Velocities(SPHERE<TV>(VECTOR<T,3>((T).5,(T).8,(T).5),(T).15),MATRIX<T,4>::Identity_Matrix(),VECTOR<T,3>());
            has_added_ball=true;}}
#endif
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
};
}
#endif
