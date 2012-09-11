//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include <climits>
#include "Armadillo_Test/ARMADILLO_TEST.h"
#include "Binding_Plasticity/BINDING_PLASTICITY_EXAMPLE.h"
#include "Binding_Springs_Test/BINDING_SPRINGS_TEST.h"
#include "Body_Test/BODY_TEST.h"
#include "Cloth_Spheres/CLOTH_SPHERES_EXAMPLE.h"
#include "Embedded_Collisions/EMBEDDED_COLLISIONS_EXAMPLE.h"
#include "Embedding_Test/EMBEDDING_TEST.h"
#include "Fragment_Tests/FRAGMENT_TESTS.h"
#include "Hair_Sim_Tests/HAIR_SIM_TESTS.h"
#include "Hair_Sim_Tests/HAIR_STRAND_TESTS.h"
#include "Hair_Tests/HAIR_TESTS.h"
#include "Incompressible/INCOMPRESSIBLE_TESTS.h"
#include "Mass_Weighted_Self_Collisions/MASS_WEIGHTED_SELF_COLLISIONS.h"
#include "Mocap_Flesh/MOCAP_FLESH.h"
#include "Rigid_Particle/RIGID_PARTICLE_EXAMPLE.h"
#include "Rigid_Particle_Mesh/RIGID_PARTICLE_MESH_EXAMPLE.h"
#include "Rigid_Particle_Spheres/RIGID_PARTICLE_SPHERES_EXAMPLE.h"
#include "Simple_Hard_Binding/SIMPLE_HARD_BINDING_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    RW rw=RW();STREAM_TYPE stream_type(rw); // gcc 3.3.2 workaround

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example;
    if(PARSE_ARGS::Find_And_Remove("-incomp",argc,argv)) example=new INCOMPRESSIBLE_TESTS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-hair_sim_tests",argc,argv)) example=new HAIR_SIM_TESTS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-hair_strand_tests",argc,argv)) example=new HAIR_STRAND_TESTS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-mass_weighted_self_collisions",argc,argv)) example=new MASS_WEIGHTED_SELF_COLLISIONS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-mocap_flesh",argc,argv)) example=new MOCAP_FLESH<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-embedded_collisions",argc,argv)) example=new EMBEDDED_COLLISIONS_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-rigid_particles",argc,argv)) example=new RIGID_PARTICLE_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-rigid_particle_mesh",argc,argv)) example=new RIGID_PARTICLE_MESH_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-rigid_particle_spheres",argc,argv)) example=new RIGID_PARTICLE_SPHERES_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-simple_hard_binding",argc,argv)) example=new SIMPLE_HARD_BINDING_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-binding_plasticity",argc,argv)) example=new BINDING_PLASTICITY_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-binding_springs",argc,argv)) example=new BINDING_SPRINGS_TEST<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-cloth_spheres",argc,argv)) example=new CLOTH_SPHERES_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-embedding_test",argc,argv)) example=new EMBEDDING_TEST<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-hair_tests",argc,argv)) example=new HAIR_TESTS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-fragment_tests",argc,argv)) example=new FRAGMENT_TESTS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-armadillo_test",argc,argv)) example=new ARMADILLO_TEST<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-body_test",argc,argv)) example=new BODY_TEST<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);

    example->want_mpi_world=true;
    PARSE_ARGS parse_args(argc,argv);
    example->Parse(parse_args);

    if(example->mpi_world->initialized) example->solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example->Adjust_Output_Directory_For_MPI(example->solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
