//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_AND_ARRAY_CONTAINER
//#####################################################################
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <Incompressible/Grid_Based_Fields/GRID_AND_ARRAY_CONTAINER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> GRID_AND_ARRAY_CONTAINER<TV,T2>::
GRID_AND_ARRAY_CONTAINER(GRID<TV>& grid_input)
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
template<class TV,class T2> GRID_AND_ARRAY_CONTAINER<TV,T2>::
~GRID_AND_ARRAY_CONTAINER()
{
    delete advection_maccormack;
    delete &advection_default;delete &boundary_default;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV,class T2> void GRID_AND_ARRAY_CONTAINER<TV,T2>::
Euler_Step(const T dt,const T time,const int number_of_ghost_cells)
{
    ARRAY<T2,TV_INT> array_ghost(grid.Cell_Indices(number_of_ghost_cells),no_init);
    boundary->Fill_Ghost_Cells(grid,array,array_ghost,dt,time,number_of_ghost_cells);
    if(cell_velocities) advection->Update_Advection_Equation_Cell(grid,array,array_ghost,*cell_velocities,*boundary,dt,time);
    else advection->Update_Advection_Equation_Cell(grid,array,array_ghost,*face_velocities,*boundary,dt,time);
}
//#####################################################################
namespace PhysBAM{
template class GRID_AND_ARRAY_CONTAINER<VECTOR<float,1>,float>;
template class GRID_AND_ARRAY_CONTAINER<VECTOR<float,2>,float>;
template class GRID_AND_ARRAY_CONTAINER<VECTOR<float,3>,float>;
template class GRID_AND_ARRAY_CONTAINER<VECTOR<double,1>,double>;
template class GRID_AND_ARRAY_CONTAINER<VECTOR<double,2>,double>;
template class GRID_AND_ARRAY_CONTAINER<VECTOR<double,3>,double>;
}
