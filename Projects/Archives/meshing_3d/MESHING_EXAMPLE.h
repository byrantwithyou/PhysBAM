//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Mike Rodgers, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MESHING_EXAMPLE
//#####################################################################
#ifndef __MESHING_EXAMPLE__
#define __MESHING_EXAMPLE__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Dyadic/READ_WRITE_DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Meshing/TETRAHEDRAL_MESHING.h>
#include "MULTIPLE_LEVELSET_IMPLICIT_SURFACE.h"
namespace PhysBAM{

template<class T>
class MESHING_EXAMPLE
{
    typedef VECTOR<T,3> TV;
public:
    const STREAM_TYPE stream_type;
    TETRAHEDRAL_MESHING<T> tetrahedral_meshing;
    std::string implicit_surface_filename;
    bool use_multiple_levelset_implicit_surface;
    bool use_octree_implicit_surface; // otherwise use level set implicit surface
    T bcc_lattice_cell_size; // size of a grid cell in the bcc lattice - 0 means it is calculated automatically
    bool use_adaptive_refinement; // otherwise uniform refinement is used
    int max_subdivision_levels;
    bool use_optimization,use_dynamics;
    int number_of_initial_optimization_steps,number_of_final_optimization_steps;
    T force_attraction_coefficient,velocity_attraction_coefficient; // for level set forces and velocities
    bool allow_tangential_velocity_slip;
    bool use_finite_element_forces; // otherwise mass/spring forces
    T youngs_modulus,poissons_ratio,rayleigh_coefficient;
    T edge_spring_stiffness,edge_spring_overdamping_fraction;
    T altitude_spring_stiffness,altitude_spring_overdamping_fraction;
    T time_step; // size of output steps - smaller time steps may be taken in the dynamics 
    int number_of_force_steps,number_of_velocity_steps;
    bool use_global_quality_criteria_for_early_exit;
    bool replace_green_refinement_with_embedded_t_junctions;
    bool allow_coarsening_to_non_graded_mesh;
    bool use_aggressive_tet_pruning;
protected:
    GRID<TV> levelset_grid; // for levelset_implicit_surface
    ARRAY<T,VECTOR<int,3> > phi3d; // for levelset_implicit_surface
    LEVELSET_IMPLICIT_OBJECT<TV> levelset_implicit_surface;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    OCTREE_GRID<T> octree_grid; // for octree_implicit_surface
#endif
    ARRAY<T> phi1d; // for octree_implicit_surface
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    DYADIC_IMPLICIT_OBJECT<TV> octree_implicit_surface;
#endif
    MULTIPLE_LEVELSET_IMPLICIT_SURFACE<T> multiple_levelset_implicit_surface;
public:

    MESHING_EXAMPLE(const STREAM_TYPE stream_type)
        :stream_type(stream_type),tetrahedral_meshing(stream_type),use_multiple_levelset_implicit_surface(false),use_octree_implicit_surface(false),
        bcc_lattice_cell_size(0),use_adaptive_refinement(true),max_subdivision_levels(7),
        use_optimization(true),use_dynamics(true),
        number_of_initial_optimization_steps(3),number_of_final_optimization_steps(0),
        force_attraction_coefficient((T).1),velocity_attraction_coefficient((T).1),allow_tangential_velocity_slip(false),
        use_finite_element_forces(false),youngs_modulus(500),poissons_ratio((T).475),rayleigh_coefficient((T).1),
        edge_spring_stiffness((T)1e-4),edge_spring_overdamping_fraction(5),altitude_spring_stiffness((T)1e-4),altitude_spring_overdamping_fraction(5),
        time_step((T).25),number_of_force_steps(40),number_of_velocity_steps(360),use_global_quality_criteria_for_early_exit(true),
        replace_green_refinement_with_embedded_t_junctions(false),allow_coarsening_to_non_graded_mesh(false),
        levelset_implicit_surface(levelset_grid,phi3d)
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        ,octree_grid(),octree_implicit_surface(octree_grid,phi1d)
#endif
    {
        tetrahedral_meshing.Use_Dynamic_Ether_Viscosity((T).5);
    }

    virtual ~MESHING_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Implicit_Surface
//#####################################################################
virtual void Initialize_Implicit_Surface()
{
    // read octree or levelset implicit surface from file
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(use_octree_implicit_surface){
        FILE_UTILITIES::Read_From_File(stream_type,implicit_surface_filename,octree_implicit_surface);
        tetrahedral_meshing.Initialize(&octree_implicit_surface);
        octree_implicit_surface.levelset.grid.Update_Maximum_Depth();
        octree_implicit_surface.Update_Minimum_Cell_Size(octree_implicit_surface.levelset.grid.maximum_depth);
        octree_implicit_surface.levelset.grid.Node_Iterator_Data();}  
    else
#endif
        if(use_multiple_levelset_implicit_surface){
        LOG::Time("reading secondary levelset files");
        FILE_UTILITIES::Create_From_File(stream_type,implicit_surface_filename,multiple_levelset_implicit_surface.primary_levelset);
        LOG::Stop_Time();
        T max_secondary_cell_size=0;
        for(int secondary_index=1;;secondary_index++){
            if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s_secondary_%d",implicit_surface_filename.c_str(),secondary_index))) break;
            multiple_levelset_implicit_surface.secondary_levelsets.Resize(secondary_index);
            LOG::Time(STRING_UTILITIES::string_sprintf("reading secondary levelset file %d",secondary_index));
            FILE_UTILITIES::Create_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s_secondary_%d",implicit_surface_filename.c_str(),secondary_index),
                multiple_levelset_implicit_surface.secondary_levelsets(secondary_index));
            LOG::Stop_Time();
            max_secondary_cell_size=max(max_secondary_cell_size,multiple_levelset_implicit_surface.secondary_levelsets(secondary_index)->levelset.grid.DX().Max());}
        tetrahedral_meshing.Initialize(&multiple_levelset_implicit_surface);
        multiple_levelset_implicit_surface.Compute_Normals();multiple_levelset_implicit_surface.Update_Box();multiple_levelset_implicit_surface.Update_Minimum_Cell_Size();
        multiple_levelset_implicit_surface.Set_Inside_Threshold((T)3*max_secondary_cell_size);}
    else{ // level set implicit surface
        FILE_UTILITIES::Read_From_File(stream_type,implicit_surface_filename,levelset_implicit_surface);
        tetrahedral_meshing.Initialize(&levelset_implicit_surface);
        levelset_implicit_surface.Compute_Normals(); // can precompute normals in this case 
        levelset_implicit_surface.Update_Box();}
}
//#####################################################################
};
}
#endif
