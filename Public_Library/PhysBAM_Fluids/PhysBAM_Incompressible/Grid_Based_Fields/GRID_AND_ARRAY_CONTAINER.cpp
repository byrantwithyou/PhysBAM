//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_AND_ARRAY_CONTAINER
//#####################################################################
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/GRID_AND_ARRAY_CONTAINER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> GRID_AND_ARRAY_CONTAINER<T_GRID,T2>::
GRID_AND_ARRAY_CONTAINER(T_GRID& grid_input)
    :grid(grid_input),advection_maccormack(0),advection_default(*new T_ADVECTION_SEMI_LAGRANGIAN_SCALAR),boundary_default(*new BOUNDARY_REFLECTION_UNIFORM<TV,T>),face_velocities(0)
{
    advection=&advection_default;
    boundary=&boundary_default;
    Initialize_Array();
    Initialize_Domain_Boundary_Conditions();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> GRID_AND_ARRAY_CONTAINER<T_GRID,T2>::
~GRID_AND_ARRAY_CONTAINER()
{
    delete advection_maccormack;
    delete &advection_default;delete &boundary_default;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID,class T2> void GRID_AND_ARRAY_CONTAINER<T_GRID,T2>::
Euler_Step(const T dt,const T time,const int number_of_ghost_cells)
{
    ARRAY<T2,TV_INT> array_ghost(grid.Cell_Indices(number_of_ghost_cells),false);boundary->Fill_Ghost_Cells(grid,array,array_ghost,dt,time,number_of_ghost_cells);
    if(cell_velocities) advection->Update_Advection_Equation_Cell(grid,array,array_ghost,*cell_velocities,*boundary,dt,time);
    else advection->Update_Advection_Equation_Cell(grid,array,array_ghost,*face_velocities,*boundary,dt,time);
}
//#####################################################################
namespace PhysBAM{
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<float,1> >,float>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<float,2> >,float>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<float,3> >,float>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<double,1> >,double>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<double,2> >,double>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<double,3> >,double>;
}
