//#####################################################################
// Copyright 2005-2017, Eran Guendelman, Craig Schroeder, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_ITERATOR
//#####################################################################
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_ITERATOR<TV>::
CELL_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells,
    const T_REGION& region_type,const int side)
    :grid(grid_input)
{
    assert(-1<=side && side<2*TV::m);
    // these types may not really make sense here    
    assert(region_type!=GRID<TV>::BOUNDARY_REGION);
    assert(region_type!=GRID<TV>::INTERIOR_REGION);

    int inner_ghost=region_type==GRID<TV>::BOUNDARY_INTERIOR_REGION?-1:0;
    this->Set_Range(grid_input.numbers_of_cells,number_of_ghost_cells,inner_ghost);
    RI flags=region_type==GRID<TV>::GHOST_REGION?RI::ghost:RI::none;
    this->Initialize(flags,side);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_ITERATOR<TV>::
CELL_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input)
    :RANGE_ITERATOR<TV::m>(region_input),grid(grid_input)
{
}
//#####################################################################
namespace PhysBAM{
template class CELL_ITERATOR<VECTOR<float,0> >;
template class CELL_ITERATOR<VECTOR<float,1> >;
template class CELL_ITERATOR<VECTOR<float,2> >;
template class CELL_ITERATOR<VECTOR<float,3> >;
template class CELL_ITERATOR<VECTOR<double,0> >;
template class CELL_ITERATOR<VECTOR<double,1> >;
template class CELL_ITERATOR<VECTOR<double,2> >;
template class CELL_ITERATOR<VECTOR<double,3> >;
}
