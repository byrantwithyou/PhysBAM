//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EDGE_ITERATOR
//##################################################################### 
#include <Grid_Tools/Grids/EDGE_ITERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
// axis_input==0 means iterate through faces in all dimensions
template<class TV> EDGE_ITERATOR<TV>::
EDGE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells_input,const T_REGION& region_type_input,const int side_input,int axis_input)
    :GRID_ITERATOR_BASE<TV>(grid_input),region_type(region_type_input),side(side_input),number_of_ghost_cells(number_of_ghost_cells_input)
{
    assert(side==-1);
    assert(region_type==GRID<TV>::WHOLE_REGION);
    if(axis_input>=0){single_axis=true;Reset_Axis(axis_input);}
    else{single_axis=false;Reset_Axis(0);}
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EDGE_ITERATOR<TV>::
EDGE_ITERATOR(const GRID<TV>& grid_input,const int axis_input,const TV_INT& face_index)
    :GRID_ITERATOR_BASE<TV>(grid_input,RANGE<TV_INT>(face_index,face_index+1))
{
    single_axis=true;
    axis=axis_input;
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EDGE_ITERATOR<TV>::
EDGE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,const int axis_input)
    :GRID_ITERATOR_BASE<TV>(grid_input),axis(axis_input),single_axis(true)
{
    assert(axis>=0);
    Add_Region(explicit_region_input);
    Reset();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void EDGE_ITERATOR<TV>::
Next_Helper()
{
    GRID_ITERATOR_BASE<TV>::Next_Helper();
    if(!valid && !single_axis && axis<TV::m-1) Reset_Axis(axis+1);
}
//#####################################################################
// Function Reset_Axis
//#####################################################################
template<class TV> void EDGE_ITERATOR<TV>::
Reset_Axis(const int axis_input)
{
    axis=axis_input;Reset_Regions();
    RANGE<TV_INT> domain(grid.Node_Indices(number_of_ghost_cells));
    switch(region_type){
        case GRID<TV>::WHOLE_REGION:
            assert(side<0);
            domain.max_corner(axis)--;
            Add_Region(domain);
            break;
        default:PHYSBAM_NOT_IMPLEMENTED();}
    Reset();
}
//#####################################################################
namespace PhysBAM{
template class EDGE_ITERATOR<VECTOR<float,1> >;
template class EDGE_ITERATOR<VECTOR<float,2> >;
template class EDGE_ITERATOR<VECTOR<float,3> >;
template class EDGE_ITERATOR<VECTOR<double,1> >;
template class EDGE_ITERATOR<VECTOR<double,2> >;
template class EDGE_ITERATOR<VECTOR<double,3> >;
}
