//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "BPP_FLUIDS_EXAMPLE.h"
#include "BPP_SOLIDS_EXAMPLE.h"
#include "FLUID_SOLVER_PB.h"
#include "PARTITIONED_DRIVER.h"
#include "SOLID_FLUID_INTERFACE_PB.h"
#include "SOLID_SOLVER_PB.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    int current_argv=0;
    ARRAY<char*> argvs[3];
    for(auto& a:argvs) a.Append(argv[0]);

    for(int i=1;i<argc;i++)
    {
        if(argv[i][0]=='-' && argv[i][1]=='-' && argv[i][2]==0) current_argv++;
        else argvs[current_argv].Append(argv[i]);
    }
    for(auto& a:argvs) a.Append(0);
    PHYSBAM_ASSERT(current_argv==2);

    PARSE_ARGS parse_args(argvs[0].m-1,argvs[0].base_pointer);
    PARSE_ARGS parse_args_s(argvs[1].m-1,argvs[1].base_pointer);
    PARSE_ARGS parse_args_f(argvs[2].m-1,argvs[2].base_pointer);

    SOLID_SOLVER_PB<TV> solid_solver;

    // SOLIDS
    {
        auto example=new BPP_SOLIDS_EXAMPLE<T>(stream_type,parse_args_s);
        example->mpi_world=new MPI_WORLD(parse_args_s);
        if(example->mpi_world->initialized) example->solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
        example->Adjust_Output_Directory_For_MPI(example->solid_body_collection.deformable_body_collection.mpi_solids);
        solid_solver.driver=new SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>(*example);
    }

    FLUID_SOLVER_PB<TV> fluid_solver;

    // FLUIDS
    {
        auto example=new BPP_FLUIDS_EXAMPLE<T>(stream_type,parse_args_f);
        example->mpi_world=new MPI_WORLD(parse_args_f);

        FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example->fluids_parameters;
        if(example->mpi_world->initialized)
            fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<TV>(*fluids_parameters.grid,3,false,VECTOR<int,2>(),fluids_parameters.periodic);
        else if(fluids_parameters.periodic!=VECTOR<bool,2>()){LOG::cerr<<"Periodic domains require MPI."<<std::endl;exit(1);}

        example->Adjust_Output_Directory_For_MPI(fluids_parameters.mpi_grid);

        fluid_solver.driver=new SOLIDS_FLUIDS_DRIVER_UNIFORM<TV>(*example);
    }

    SOLID_FLUID_INTERFACE_PB<TV> interface;

    PARTITIONED_DRIVER<TV> dr;
    dr.solid_solver=&solid_solver;
    dr.fluid_solver=&fluid_solver;
    dr.interface=&interface;

    parse_args.Add("-last_frame",&dr.last_frame,"frame","number of frames to simulate");
    parse_args.Add("-frame_dt",&dr.frame_dt,"dt","Frame size");
    parse_args.Add("-min_dt",&dr.min_dt,"dt","Minimum time step size");
    parse_args.Add("-max_dt",&dr.max_dt,"dt","Maximum time step size");
    parse_args.Add("-dt",&dr.fixed_dt,"dt","Maximum time step size");
    parse_args.Add("-max_subit",&dr.max_subiterations,"num","Maximum time step size");
    parse_args.Add("-utol",&dr.utol,"tol","Fluid velocity tolerance");
    parse_args.Add("-ptol",&dr.ptol,"tol","Fluid pressure tolerance");
    parse_args.Add("-xtol",&dr.xtol,"tol","Solids position tolerance");
    parse_args.Add("-vtol",&dr.vtol,"tol","Solids velocity tolerance");
    parse_args.Add("-p0tol",&dr.p0tol,"tol","BPP p0 tolerance");
    parse_args.Parse();
    dr.Run();

    delete &solid_solver.driver->example;
    delete solid_solver.driver;
    delete &fluid_solver.driver->example;
    delete fluid_solver.driver;

    return 0;
}
//#####################################################################
