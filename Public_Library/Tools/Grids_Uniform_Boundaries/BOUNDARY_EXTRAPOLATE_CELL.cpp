//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EXTRAPOLATE_CELL
//#####################################################################
#include <Tools/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_EXTRAPOLATE_CELL.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_EXTRAPOLATE_CELL<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost); // interior
    RANGE<TV_INT> domain=grid.Domain_Indices(number_of_ghost_cells);
    EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::Extrapolate_Cell(grid,
        [&](const TV_INT& index){return domain.Lazy_Inside_Half_Open(index);},
        number_of_ghost_cells,u_ghost,3,number_of_ghost_cells);
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY_EXTRAPOLATE_CELL<TV,T2>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const
{
}
namespace PhysBAM{
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<float,1>,float>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<float,2>,float>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<float,3>,float>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<double,1>,double>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<double,2>,double>;
template class BOUNDARY_EXTRAPOLATE_CELL<VECTOR<double,3>,double>;
}
