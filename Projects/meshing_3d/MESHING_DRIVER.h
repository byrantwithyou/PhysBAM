//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Mike Rodgers, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MESHING_DRIVER
//#####################################################################
#ifndef __MESHING_DRIVER__
#define __MESHING_DRIVER__

#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Dynamics/Meshing/TETRAHEDRAL_MESHING.h>
#include "MESHING_EXAMPLE.h"
using namespace PhysBAM;

template<class T>
class MESHING_DRIVER
{
    typedef VECTOR<T,3> TV;
public:
    MESHING_EXAMPLE<T>& example;

    MESHING_DRIVER(MESHING_EXAMPLE<T>& example_input)
        :example(example_input)
    {
        Initialize();
    }

//#####################################################################
// Function Initialize
//#####################################################################
void Initialize()
{
    LOG::SCOPE scope("INITIALIZATION","initializing");
    // implicit surface
    example.Initialize_Implicit_Surface();

    // tetrahedral meshing
    example.tetrahedral_meshing.Replace_Green_Refinement_With_Embedded_T_Junctions(example.replace_green_refinement_with_embedded_t_junctions,example.allow_coarsening_to_non_graded_mesh);
    if(example.use_dynamics){
        example.tetrahedral_meshing.level_set_forces_and_velocities->Set_Force_Attraction_Coefficient(example.force_attraction_coefficient);
        example.tetrahedral_meshing.level_set_forces_and_velocities->Set_Velocity_Attraction_Coefficient(example.velocity_attraction_coefficient);
        example.tetrahedral_meshing.level_set_forces_and_velocities->Allow_Tangential_Velocity_Slip(example.allow_tangential_velocity_slip);
        if(example.use_finite_element_forces) example.tetrahedral_meshing.Use_Finite_Elements(example.youngs_modulus,example.poissons_ratio,example.rayleigh_coefficient);
        else example.tetrahedral_meshing.Use_Masses_And_Springs(example.edge_spring_stiffness,example.edge_spring_overdamping_fraction,
            example.altitude_spring_stiffness,example.altitude_spring_overdamping_fraction);}
    if(example.use_optimization)
        example.tetrahedral_meshing.Use_Global_Quality_Criteria_For_Early_Exit(example.use_global_quality_criteria_for_early_exit);
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
void Execute_Main_Program()
{  
    example.tetrahedral_meshing.Create_Initial_Mesh(example.bcc_lattice_cell_size,example.use_adaptive_refinement,example.max_subdivision_levels,true,false,example.use_aggressive_tet_pruning);
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=example.tetrahedral_meshing.solid_body_collection.deformable_body_collection;
    LOG::cout<<"tetrahedrons = "<<deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>().mesh.elements.m<<std::endl;
    LOG::cout<<"particles = "<<deformable_body_collection.particles.Size()<<std::endl;

    if(example.use_optimization){
        example.tetrahedral_meshing.Initialize_Optimization();
        example.tetrahedral_meshing.Create_Final_Mesh_With_Optimization(example.number_of_initial_optimization_steps,example.number_of_final_optimization_steps);}
    if(example.use_dynamics){
        example.tetrahedral_meshing.Initialize_Dynamics();
        example.tetrahedral_meshing.Create_Final_Mesh_With_Dynamics(example.time_step,example.number_of_force_steps,example.number_of_velocity_steps);}
    example.tetrahedral_meshing.Snap_Nodes_To_Level_Set_Boundary();
    FILE_UTILITIES::Write_To_File(example.stream_type,example.tetrahedral_meshing.output_directory,deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>());
}
//##############################################################################
};
#endif

