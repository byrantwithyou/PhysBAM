//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ITERATOR
//##################################################################### 
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
// axis_input==0 means iterate through faces in all dimensions
template<class TV> FACE_ITERATOR<TV>::
FACE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells_input,
    const typename GRID<TV>::REGION& region_type_input,const int side_input,
    int axis_input)
    :grid(grid_input)
{
    int inner_ghost=0;
    if(region_type_input==GRID<TV>::BOUNDARY_REGION) inner_ghost=number_of_ghost_cells_input;
    Set_Range(grid_input.numbers_of_cells,number_of_ghost_cells_input,inner_ghost);
    RF flags=RF::none;
    switch(region_type_input){
        case GRID<TV>::WHOLE_REGION:flags=RF::interior;break;
        case GRID<TV>::GHOST_REGION:flags=RF::skip_inner;break;
        case GRID<TV>::BOUNDARY_REGION:flags=RF::none;break;
        case GRID<TV>::INTERIOR_REGION:flags=RF::interior|RF::skip_outer;break;
        case GRID<TV>::BOUNDARY_INTERIOR_REGION:PHYSBAM_FATAL_ERROR();break;}
    Initialize(flags,side_input,axis_input);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FACE_ITERATOR<TV>::
FACE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,
    const int axis_input)
    :FACE_RANGE_ITERATOR<TV::m>(explicit_region_input,RF::interior,axis_input),
    grid(grid_input)
{
}
//#####################################################################
namespace PhysBAM{
template class FACE_ITERATOR<VECTOR<float,1> >;
template class FACE_ITERATOR<VECTOR<float,2> >;
template class FACE_ITERATOR<VECTOR<float,3> >;
template class FACE_ITERATOR<VECTOR<double,1> >;
template class FACE_ITERATOR<VECTOR<double,2> >;
template class FACE_ITERATOR<VECTOR<double,3> >;
}
