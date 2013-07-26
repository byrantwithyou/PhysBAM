//#####################################################################
// Copyright 2002-2009, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Incompressible/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_FACE.h>
#include <Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Incompressible_Fluid_Evolution/INCOMPRESSIBLE_FLUID_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
INCOMPRESSIBLE_FLUID_EVOLUTION(const GRID<TV>& grid_input)
    :grid(grid_input.Get_MAC_Grid()),max_time_step(100),boundary_default(*new BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>)
{ 
    boundary=&boundary_default;
    Initialize_Grids(grid);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
~INCOMPRESSIBLE_FLUID_EVOLUTION()
{
    delete &boundary_default;
    fluids_forces.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Advance_One_Time_Step_Convection
//#####################################################################
template<class TV> void INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
Advance_One_Time_Step_Convection(const T dt,const T time,const T_FACE_ARRAYS_SCALAR& advecting_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect,const int number_of_ghost_cells)
{
    // TODO: make efficient if advection velocities are same as advected velocities
    // find ghost cells
    T_FACE_ARRAYS_SCALAR advection_face_velocities_ghost(grid,number_of_ghost_cells,false);
    T_FACE_ARRAYS_SCALAR face_velocities_to_advect_ghost(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,time,number_of_ghost_cells);
    boundary->Fill_Ghost_Faces(grid,advecting_face_velocities,advection_face_velocities_ghost,time,number_of_ghost_cells);
    
    // update convection
    advection->Update_Advection_Equation_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,advection_face_velocities_ghost,*boundary,dt,time);
}
//#####################################################################
// Function Advance_One_Time_Step_Forces
//#####################################################################
template<class TV> void INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const int number_of_ghost_cells)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);
    for(int k=0;k<fluids_forces.m;k++) fluids_forces(k)->Add_Explicit_Forces(grid,face_velocities_ghost,face_velocities,dt,time);
    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class TV> void INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
Advance_One_Time_Step_Implicit_Part(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    // boundary conditions
    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
    
    // this check moves inside incompressibility force
    //assert(Consistent_Boundary_Conditions(face_velocities));

    int ghost_cells=3;
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(grid,ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time,ghost_cells);
    for(int k=0;k<fluids_forces.m;k++) fluids_forces(k)->Add_Implicit_Forces_Before_Projection(grid,face_velocities_ghost,face_velocities,dt,time);
    for(int k=0;k<fluids_forces.m;k++) fluids_forces(k)->Add_Implicit_Forces_Projection(grid,face_velocities_ghost,face_velocities,dt,time);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
CFL(T_FACE_ARRAYS_SCALAR& face_velocities) const
{
    // advection should define its own cfl?
    T dt_convection=0;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        T local_V_norm=0;
        for(int axis=0;axis<TV::m;axis++)
            local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),
                face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}
    T dt_force=0;
    for(int k=0;k<fluids_forces.m;k++) dt_force+=fluids_forces(k)->CFL(grid,face_velocities);
    T dt_overall=(dt_convection+sqrt(sqr(dt_convection)+4*dt_force))/2;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class TV> void INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
Initialize_Grids(const GRID<TV>& grid_input)
{
    grid=grid_input.Get_MAC_Grid();
    for(int k=0;k<fluids_forces.m;k++) fluids_forces(k)->Initialize_Grids(grid);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int INCOMPRESSIBLE_FLUID_EVOLUTION<TV>::
Add_Force(INCOMPRESSIBLE_FLUIDS_FORCES<TV>* force)
{
    fluids_forces.Append(force);
    return fluids_forces.Size();
}
//#####################################################################
namespace PhysBAM{
template class INCOMPRESSIBLE_FLUID_EVOLUTION<VECTOR<float,1> >;
template class INCOMPRESSIBLE_FLUID_EVOLUTION<VECTOR<float,2> >;
template class INCOMPRESSIBLE_FLUID_EVOLUTION<VECTOR<float,3> >;
template class INCOMPRESSIBLE_FLUID_EVOLUTION<VECTOR<double,1> >;
template class INCOMPRESSIBLE_FLUID_EVOLUTION<VECTOR<double,2> >;
template class INCOMPRESSIBLE_FLUID_EVOLUTION<VECTOR<double,3> >;
}
