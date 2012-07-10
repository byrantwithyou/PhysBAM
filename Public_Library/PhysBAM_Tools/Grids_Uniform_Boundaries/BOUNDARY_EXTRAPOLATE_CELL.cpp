//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EXTRAPOLATE_CELL
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_EXTRAPOLATE_CELL.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_EXTRAPOLATE_CELL<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,TV::m> >& u,ARRAYS_ND_BASE<VECTOR<T2,TV::m> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    ARRAYS_ND_BASE<VECTOR<T2,TV::m> >::Put(u,u_ghost); // interior
    ARRAY<bool,TV_INT> inside_mask(grid.Domain_Indices(number_of_ghost_cells));
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        inside_mask(it.index)=true;

    EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::Extrapolate_Cell(grid,inside_mask,number_of_ghost_cells,u_ghost,3,number_of_ghost_cells);
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY_EXTRAPOLATE_CELL<TV,T2>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<VECTOR<T2,TV::m> >& u,const T time)
{
}
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<float,1>,float>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<float,2>,float>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<float,3>,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<double,1>,double>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<double,2>,double>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<double,3>,double>;
#endif
