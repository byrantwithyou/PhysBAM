//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_INTERIOR_ITERATOR_CELL
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_INTERIOR_ITERATOR_CELL.h>
using namespace PhysBAM;
//#####################################################################
// Function UNIFORM_INTERIOR_ITERATOR_CELL
//#####################################################################
template<class TV> UNIFORM_INTERIOR_ITERATOR_CELL<TV>::
UNIFORM_INTERIOR_ITERATOR_CELL(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input,
    const typename GRID<TV>::REGION& region_type_input,const int side_input)
    :BASE(info.grid,number_of_ghost_cells_input,region_type_input,side_input),outside_fluid(*info.outside_fluid)
{
}
//#####################################################################
namespace PhysBAM{
template class UNIFORM_INTERIOR_ITERATOR_CELL<VECTOR<float,1> >;
template class UNIFORM_INTERIOR_ITERATOR_CELL<VECTOR<float,2> >;
template class UNIFORM_INTERIOR_ITERATOR_CELL<VECTOR<float,3> >;
template class UNIFORM_INTERIOR_ITERATOR_CELL<VECTOR<double,1> >;
template class UNIFORM_INTERIOR_ITERATOR_CELL<VECTOR<double,2> >;
template class UNIFORM_INTERIOR_ITERATOR_CELL<VECTOR<double,3> >;
}
