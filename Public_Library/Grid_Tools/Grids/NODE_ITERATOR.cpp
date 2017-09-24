//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NODE_ITERATOR
//#####################################################################
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Function NODE_ITERATOR
//#####################################################################
template<class TV> NODE_ITERATOR<TV>::
NODE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells,
    const T_REGION& region_type,const int side)
    :grid(grid_input)
{
    assert(-1<=side&&side<2*TV::m);
    assert(region_type!=GRID<TV>::BOUNDARY_INTERIOR_REGION);

    int inner_ghost=0,outer_ghost=number_of_ghost_cells;
    RI flags=RI::none;
    switch(region_type){
        case GRID<TV>::WHOLE_REGION:break;
        case GRID<TV>::GHOST_REGION:flags=RI::ghost;break;
        case GRID<TV>::BOUNDARY_REGION:inner_ghost=outer_ghost-1;flags=RI::ghost;break;
        case GRID<TV>::INTERIOR_REGION:outer_ghost--;break;
        case GRID<TV>::BOUNDARY_INTERIOR_REGION:PHYSBAM_FATAL_ERROR();break;}
    this->Set_Range(grid_input.Numbers_Of_Nodes(),outer_ghost,inner_ghost);
    this->Initialize(flags,side);
}
//#####################################################################
// Function NODE_ITERATOR
//#####################################################################
template<class TV> NODE_ITERATOR<TV>::
NODE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input)
    :RANGE_ITERATOR<TV::m>(region_input),grid(grid_input)
{
}
//#####################################################################
namespace PhysBAM{
template class NODE_ITERATOR<VECTOR<float,1> >;
template class NODE_ITERATOR<VECTOR<float,2> >;
template class NODE_ITERATOR<VECTOR<float,3> >;
template class NODE_ITERATOR<VECTOR<double,1> >;
template class NODE_ITERATOR<VECTOR<double,2> >;
template class NODE_ITERATOR<VECTOR<double,3> >;
}
