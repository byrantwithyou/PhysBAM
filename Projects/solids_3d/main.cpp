//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
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
#include "Rigid_Particle/RIGID_PARTICLE_EXAMPLE.h"
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

    bool opt_incomp=false,opt_hair_sim_tests=false,opt_hair_strand_tests=false,opt_mass_weighted_self_collisions=false;
    bool opt_embedded_collisions=false,opt_rigid_particles=false;
    bool opt_rigid_particle_spheres=false,opt_simple_hard_binding=false,opt_binding_plasticity=false,opt_binding_springs=false;
    bool opt_cloth_spheres=false,opt_embedding_test=false,opt_hair_tests=false,opt_fragment_tests=false,opt_armadillo_test=false;
    bool opt_body_test=false;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-incomp",&opt_incomp,"Use incomp test");
    parse_args.Add("-hair_sim_tests",&opt_hair_sim_tests,"Use hair_sim_tests test");
    parse_args.Add("-hair_strand_tests",&opt_hair_strand_tests,"Use hair_strand_tests test");
    parse_args.Add("-mass_weighted_self_collisions",&opt_mass_weighted_self_collisions,"Use mass_weighted_self_collisions test");
    parse_args.Add("-embedded_collisions",&opt_embedded_collisions,"Use embedded_collisions test");
    parse_args.Add("-rigid_particles",&opt_rigid_particles,"Use rigid_particles test");
    parse_args.Add("-rigid_particle_spheres",&opt_rigid_particle_spheres,"Use rigid_particle_spheres test");
    parse_args.Add("-simple_hard_binding",&opt_simple_hard_binding,"Use simple_hard_binding test");
    parse_args.Add("-binding_plasticity",&opt_binding_plasticity,"Use binding_plasticity test");
    parse_args.Add("-binding_springs",&opt_binding_springs,"Use binding_springs test");
    parse_args.Add("-cloth_spheres",&opt_cloth_spheres,"Use cloth_spheres test");
    parse_args.Add("-embedding_test",&opt_embedding_test,"Use embedding_test test");
    parse_args.Add("-hair_tests",&opt_hair_tests,"Use hair_tests test");
    parse_args.Add("-fragment_tests",&opt_fragment_tests,"Use fragment_tests test");
    parse_args.Add("-armadillo_test",&opt_armadillo_test,"Use armadillo_test test");
    parse_args.Add("-body_test",&opt_body_test,"Use body_test test");
    parse_args.Parse(true);

    SOLIDS_EXAMPLE<TV>* example;
    if(opt_incomp) example=new INCOMPRESSIBLE_TESTS<T>(stream_type);
    else if(opt_hair_sim_tests) example=new HAIR_SIM_TESTS<T>(stream_type);
    else if(opt_hair_strand_tests) example=new HAIR_STRAND_TESTS<T>(stream_type);
    else if(opt_mass_weighted_self_collisions) example=new MASS_WEIGHTED_SELF_COLLISIONS<T>(stream_type);
    else if(opt_embedded_collisions) example=new EMBEDDED_COLLISIONS_EXAMPLE<T>(stream_type);
    else if(opt_rigid_particles) example=new RIGID_PARTICLE_EXAMPLE<T>(stream_type);
    else if(opt_rigid_particle_spheres) example=new RIGID_PARTICLE_SPHERES_EXAMPLE<T>(stream_type);
    else if(opt_simple_hard_binding) example=new SIMPLE_HARD_BINDING_EXAMPLE<T>(stream_type);
    else if(opt_binding_plasticity) example=new BINDING_PLASTICITY_EXAMPLE<T>(stream_type);
    else if(opt_binding_springs) example=new BINDING_SPRINGS_TEST<T>(stream_type);
    else if(opt_cloth_spheres) example=new CLOTH_SPHERES_EXAMPLE<T>(stream_type);
    else if(opt_embedding_test) example=new EMBEDDING_TEST<T>(stream_type);
    else if(opt_hair_tests) example=new HAIR_TESTS<T>(stream_type);
    else if(opt_fragment_tests) example=new FRAGMENT_TESTS<T>(stream_type);
    else if(opt_armadillo_test) example=new ARMADILLO_TEST<T>(stream_type);
    else if(opt_body_test) example=new BODY_TEST<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);

    example->want_mpi_world=true;
    example->Parse(parse_args);

    if(example->mpi_world->initialized) example->solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example->Adjust_Output_Directory_For_MPI(example->solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
