//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Incompressible/Forces/VISCOSITY.h>
#include <Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VISCOSITY<TV>::
VISCOSITY(LAPLACE_UNIFORM<TV>& elliptic_solver_input,const ARRAY<T,TV_INT>& variable_viscosity_input,const T density_input,const T viscosity_input,bool implicit_viscosity_input,
    bool use_explicit_part_of_implicit_viscosity_input,bool use_variable_viscosity_input,int maximum_implicit_viscosity_iterations_input,bool use_psi_R_input)
    :elliptic_solver(elliptic_solver_input),density(density_input),viscosity(viscosity_input),implicit_viscosity(implicit_viscosity_input),
    use_explicit_part_of_implicit_viscosity(use_explicit_part_of_implicit_viscosity_input),variable_viscosity(variable_viscosity_input),use_variable_viscosity(use_variable_viscosity_input),
    maximum_implicit_viscosity_iterations(maximum_implicit_viscosity_iterations_input),use_psi_R(use_psi_R_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> VISCOSITY<TV>::
~VISCOSITY()
{
}
//#####################################################################
// Function Add_Explicit_Forces
//#####################################################################
template<class TV> void VISCOSITY<TV>::
Add_Explicit_Forces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    if(dt && use_variable_viscosity && (!implicit_viscosity || use_explicit_part_of_implicit_viscosity)){
        if(!implicit_viscosity) PHYSBAM_NOT_IMPLEMENTED();
        IMPLICIT_VISCOSITY_UNIFORM<TV>::Variable_Viscosity_Explicit_Part(density,variable_viscosity,grid,face_velocities,face_velocities_ghost,dt,time);}
}
//#####################################################################
// Function Add_Implicit_Forces_Before_Projection
//#####################################################################
template<class TV> void VISCOSITY<TV>::
Add_Implicit_Forces_Before_Projection(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    if(!dt || (!use_variable_viscosity && viscosity==0)) return;
    for(int axis=0;axis<TV::m;axis++){
        IMPLICIT_VISCOSITY_UNIFORM<TV> implicit_viscosity(elliptic_solver,variable_viscosity,density,viscosity,elliptic_solver.mpi_grid,axis,use_variable_viscosity,use_psi_R);
        implicit_viscosity.Viscous_Update(grid,face_velocities,face_velocities_ghost,dt,time,maximum_implicit_viscosity_iterations);}
    if(elliptic_solver.mpi_grid) elliptic_solver.mpi_grid->Copy_Common_Face_Data(face_velocities);
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class TV> void VISCOSITY<TV>::
Initialize_Grids(const GRID<TV>& grid)
{
    elliptic_solver.Initialize_Grid(grid);
}
//#####################################################################
namespace PhysBAM{
template class VISCOSITY<VECTOR<float,1> >;
template class VISCOSITY<VECTOR<float,2> >;
template class VISCOSITY<VECTOR<float,3> >;
template class VISCOSITY<VECTOR<double,1> >;
template class VISCOSITY<VECTOR<double,2> >;
template class VISCOSITY<VECTOR<double,3> >;
}
