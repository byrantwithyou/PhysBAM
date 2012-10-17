//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/GRADIENT_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Dynamics/Drivers/COMPRESSIBLE_EXAMPLE.h>
using namespace PhysBAM;
//#####################################################################
// COMPRESSIBLE_EXAMPLE
//#####################################################################
template<class TV> COMPRESSIBLE_EXAMPLE<TV>::
COMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input)
    :EXAMPLE<TV>(stream_type_input),
    number_of_ghost_cells(3),cfl((T).9),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    euler(mac_grid),euler_solid_fluid_coupling_utilities(euler),compressible_fluid_collection(mac_grid),face_velocities(mac_grid),
    conservation_method(0),boundary(0),pressure_boundary(0),
    collision_bodies_affecting_fluid(new T_GRID_BASED_COLLISION_GEOMETRY(mac_grid)),
    left_wall(true),right_wall(true),bottom_wall(true),top_wall(false),front_wall(true),back_wall(true),timesplit(false),
    set_max_time_step(false),max_time_step(1.e8),spatial_order(3),
    rungekutta_order(3),tolerance((T)1e-8),iterations(0),solve_single_cell_neumann_regions(false),
    fluid_affects_solid(false),solid_affects_fluid(false),apply_isobaric_fix(false),monitor_conservation_error(false),
    perform_rungekutta_for_implicit_part(false),use_sound_speed_for_cfl(false),
    use_sound_speed_based_dt_multiple_for_cfl(false),multiplication_factor_for_sound_speed_based_dt(1.)
{}
//#####################################################################
// ~COMPRESSIBLE_EXAMPLE
//#####################################################################
template<class TV> COMPRESSIBLE_EXAMPLE<TV>::
~COMPRESSIBLE_EXAMPLE()
{}
//#####################################################################
// Function Apply_Isobaric_Fix
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Apply_Isobaric_Fix(const T dt,const T time)
{
    if(apply_isobaric_fix) euler_solid_fluid_coupling_utilities.Apply_Isobaric_Fix(dt,time);
}
//#####################################################################
// Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Set_Domain_Boundary_Conditions()
{}
//#####################################################################
// Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    euler.euler_projection.Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Set_Neumann_Boundary_Conditions
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Set_Neumann_Boundary_Conditions()
{
    collision_bodies_affecting_fluid->Compute_Psi_N(euler.euler_projection.elliptic_solver->psi_N,&face_velocities);
}
//#####################################################################
// Set_Boundary_Conditions
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time)
{
    euler.euler_projection.elliptic_solver->psi_N.Fill(false);euler.euler_projection.elliptic_solver->psi_D.Fill(false);
    Set_Domain_Boundary_Conditions();
    Set_Dirichlet_Boundary_Conditions(time);
    Set_Neumann_Boundary_Conditions();
}
//#####################################################################
// Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Initialize_Solid_Fluid_Coupling()
{
    euler_solid_fluid_coupling_utilities.Initialize_Solid_Fluid_Coupling(collision_bodies_affecting_fluid);
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",mac_grid);
    compressible_fluid_collection.Write_Output_Files(stream_type,output_directory,frame);

    if(euler.timesplit) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/compressible_implicit_pressure",euler.euler_projection.p);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",euler.euler_projection.elliptic_solver->psi_D);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",euler.euler_projection.elliptic_solver->psi_N);
}
//#####################################################################
template class COMPRESSIBLE_EXAMPLE<VECTOR<float,1> >;
template class COMPRESSIBLE_EXAMPLE<VECTOR<float,2> >;
template class COMPRESSIBLE_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_EXAMPLE<VECTOR<double,1> >;
template class COMPRESSIBLE_EXAMPLE<VECTOR<double,2> >;
template class COMPRESSIBLE_EXAMPLE<VECTOR<double,3> >;
#endif
