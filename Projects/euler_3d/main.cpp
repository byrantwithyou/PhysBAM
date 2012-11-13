//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Incompressible_Drop/INCOMPRESSIBLE_DROP.h"
#include "Sphere_Example/SPHERE_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"

using namespace PhysBAM;

template<class T> void main_program(PARSE_ARGS& parse_args){
    typedef T RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,3> TV;

    bool incompressible=false,opt_sphere=false,opt_drop=false;
    VECTOR<int,3> procs;
    parse_args.Add("-incompressible",&incompressible,"Incompressible");
    parse_args.Add("-sphere",&opt_sphere,"Sphere");
    parse_args.Add("-drop",&opt_drop,"Drop");
    parse_args.Add("-xprocs",&procs.x,"n","X procs");
    parse_args.Add("-yprocs",&procs.y,"n","Y procs");
    parse_args.Add("-zprocs",&procs.z,"n","Z procs");
    parse_args.Parse();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(opt_sphere) example=new SPHERE_EXAMPLE<T>(stream_type);
    else if(opt_drop) example=new INCOMPRESSIBLE_DROP<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type,incompressible);
    example->want_mpi_world=true;
    example->Parse(parse_args);

    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*example->fluids_parameters.grid,3,false,procs,
                VECTOR<bool,3>(),example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
            example->solid_body_collection.deformable_body_collection.simulate=false;
            example->solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;}
        else{
            example->fluids_parameters.simulate=false;
            example->solids_fluids_parameters.mpi_solid_fluid->Create_Fluid_Comm_For_Solid_Nodes();}
    }
    example->Adjust_Output_Directory_For_MPI(example->solids_fluids_parameters.mpi_solid_fluid);
    
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();
    
    delete example;
}

int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    bool type_double=false;
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Parse(true);

    if(type_double)
        main_program<double>(parse_args);
    else
        main_program<float>(parse_args);
    return 0;
}
//#####################################################################

